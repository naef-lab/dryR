#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param countData	matrix containing non-negative integers; each column represents a sample, each row represents a gene/transcript.
#' @param group	vector containing the name of each sample.
#' @param time	vector containing numeric values of the time for each sample.
#' @param countData	matrix containing non-negative integers; each column represents a sample, each row represents a gene/transcript.
#' @param T_	numeric value to indicate period length of the oscillation.
#' @param batch	vector containing potential batch effects between samples.
#' @param nthreads	vector numeric value to indicate the threads for parallel computing .
#' @export
#' @examples load("simulatedData.RData")
#' @examples countData = simData$abundData
#' dryseq_function()
dryseq=function(countData,group,time,T_=24,sample_name=colnames(countData),batch=rep("A",length(sample_name)),n.cores=round(detectCores()*.6,0) ){
  require('DESeq2')
  require("combinat")
  require("parallel")
  require("gplots")
  library("RColorBrewer")
  library("doMC")
  library("circular")

  register(MulticoreParam(n.cores))
  registerDoMC(cores=n.cores)

  countData = countData[rowSums(countData)!=0,]

  s1 <- sin(2*pi*time/T_)
  c1 <- cos(2*pi*time/T_)

  conds  = cbind(group,s1,c1,batch)
  colnames(conds) = c("group","s1","c1","batch")

  colData <- data.frame(row.names=colnames(countData), conds)
  N=length(unique(group))

  ############################
  # FIT RHYTHMS

  message("fitting rhythmic models")

  models = create_matrix_list(time,N,T_)
  #Reorder u, a, b
  models = lapply(models, function(l) l[,c(grep("u",colnames(l)),grep("a|b",colnames(l)))]   )

  for (i in 1:length(models)){
    rownames(models[[i]]) = rownames(colData)}

  if (length(unique(batch))>1) {
    # add the batch effect
    model_b = as.matrix(model.matrix(~  batch),contrasts.arg=NULL)[,2:length(unique(batch)),drop=F]
    colnames(model_b)=paste0("BATCH_",unique(batch)[-1])
    models = lapply(models, function(l) cbind(model_b,l))
    models = lapply(models, function(l) l[,c(grep("u",colnames(l)),grep("BATCH",colnames(l)),grep("a|b",colnames(l)))]   )
  }



  dds = DESeqDataSetFromMatrix(countData = countData, colData = colData,   design = models[[length(models)]])
  dds.full = DESeq(dds, full=models[[length(models)]], betaPrior = F, fitType = "parametric", test = "Wald", parallel =T, quiet = T)

  deviances = sapply(models[-length(models)], function(m){
    dds.x = nbinomWaldTest(dds.full, modelMatrix= m, betaPrior = F, quiet = T)
    return(mcols(dds.x)$deviance)
  })

  deviances = cbind(deviances,mcols(dds.full)$deviance)

  message("computing BICW (rhythm)")

  # calculate the BIC
  BIC = as.data.frame(sapply(1:ncol(deviances), function(i) { deviances[,i] + log(ncol(countData)) * ncol(models[[i]] )}   ))

  #calculate the BICW
  # I think that the Zaehler is actually 1 because the difference between the model j (=choosen model) and the minimum model is 0
  BICW               = t(apply(BIC,1,compute_BICW))
  choosen_model      = apply(BIC,1,which.min)
  choosen_model_BICW = apply(BICW,1,max)

  ############################
  # FIT BASELINE

  message("fitting mean models")

  model_mean_cond=create_matrix_list_mean(N,group)
  model_mean_cond=lapply(model_mean_cond,annotate_matrix,group)

  for (i in 1:length(model_mean_cond)){
    rownames(model_mean_cond[[i]]) = rownames(colData)}

  # choose the best model for rhthmicity and then run the mean on the samples

  #dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = models[[length(models)]])
  #dds.full = DESeq(dds, full=models[[length(models)]], betaPrior = F, fitType = "parametric", test = "Wald")


  DDS_dev =  foreach (i = 1:length(models)) %dopar% {  #nrow(dds.full)

    sel = which(choosen_model==i)
    gene = rownames(dds.full)[sel]


    if(length(gene)>0){
      M=models[[i]]
      #build the gene specific model from the rhythmic point of view
      gene_specific_mean_models = lapply(model_mean_cond,
                                         function(x) cbind(x,M[,-grep("u",colnames(M))]))

      dev <- lapply(gene_specific_mean_models,function(m){
        dds.m <- dds.full # Copying the full model
        dds.m <- nbinomWaldTest(dds.m[gene], modelMatrix= as.matrix(m), betaPrior = F) # Re-run wald test
        #return(list(dds.m,mcols(dds.m)$deviance)) # Returning deviances (-2 * log likelihood) // https://support.bioconductor.org/p/107472/
        return(list(dds.m, mcols(dds.m)$deviance)) # Returning deviances (-2 * log likelihood) // https://support.bioconductor.org/p/107472/

      })
    }

    if(length(gene)==0){dev = list (NA, NA)}



    return(dev)
  }

  #DDS_dev: 1st list: cm rhythmicicty |2nd list cm_mean | 3rd list dds or deviance ||
  #DDS_dev[[3]][[3]][[2]]
  #rownames(DDS_dev[[1]][[1]][[1]])

  deviance_mean = NULL
  for (cm_r in 1:length(models)){

    if(!is.na(DDS_dev[[cm_r]][1])){
      deviance_mean.x  = rbind(sapply(1:15,function(x) {DDS_dev[[cm_r]][[x]][[2]]}))
      rownames(deviance_mean.x)  = rownames(DDS_dev[[cm_r]][[1]][[1]])
      deviance_mean    = rbind(deviance_mean, deviance_mean.x)}
  }


  # ordnen nach counts_Data reihenfolge!!!
  deviance_mean = deviance_mean[rownames(countData),]

  message("computing BICW (mean)")


  # calculate the BIC
  BIC_mean = as.data.frame(sapply(1:ncol(deviance_mean), function(i) { deviance_mean[,i] + log(ncol(countData)) * ncol(model_mean_cond[[i]] )}   ))

  #calculate the BICW
  BICW_mean = t(apply(BIC_mean,1,compute_BICW))

  choosen_model_mean = apply(BIC_mean,1,which.min)
  choosen_model_mean_BICW = apply(BICW_mean,1,max)

  ################
  # coefficients / mean, amplitude and phase
  ###############

  message("extracting rhythmic parameters")

  parameters=NULL

  parameters =  foreach (i = 1:nrow(deviance_mean)) %dopar% {
    gene = rownames(deviance_mean)[i]
    cm_r = choosen_model[i]
    cm_m = choosen_model_mean[i]
    dds= DDS_dev[[cm_r]][[cm_m]][[1]]
    out = compute_param(dds, gene ,T_,N)
    return(data.frame(row.names= gene, t(matrix(out)))           )
  }

  parameters            = data.frame(do.call(rbind.data.frame, parameters))
  colnames(parameters)  = c(paste(c('mean','a','b','amp','relamp','phase'),rep(unique(group),each =6), sep = "_"))
  parameters            = parameters[rownames(countData),]

  # Generate all the count and expression data
  # raw counts
  counts_RF        =  counts(dds.full, normalized = FALSE)

  #vst stabilized counts
  vsd <- varianceStabilizingTransformation(dds.full)
  vsd <- assay(vsd)

  #normalized counts
  ncounts_RF       = counts(dds.full, normalized = TRUE)

  # generate a table summarizing the analysis

  complete_parametes = cbind(parameters,choosen_model,choosen_model_BICW, choosen_model_mean, choosen_model_mean_BICW)
  global_table = merge(ncounts_RF,complete_parametes, by="row.names")
  rownames(global_table) = global_table$Row.names
  global_table_df  = global_table[,-grep("Row.names",colnames(global_table))]

  global_table_df = global_table_df[rownames(countData),]

  out = list()

  out[["results"]]     = global_table_df
  out[["BICW_rhythm"]] = BICW
  out[["BICW_mean"]]   = BICW_mean
  out[["vsd"]]         = vsd
  out[["ncounts"]]     = ncounts_RF
  out[["counts"]]      = counts_RF
  out[["parameters"]]  = complete_parametes
  out[["cook"]]        = assays(dds.full)[["cooks"]]
  out[["dds.full"]]    = dds.full

  message("finished!")
  return(out)
  # to add flags for low expression, high cook's distance

}
