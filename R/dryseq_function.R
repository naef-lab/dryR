#' Differential rhythmicity analysis for RNA-Seq datasets
#'
#' This function performs a rhythmicity analysis based on generalized linear models with a subsequent models selection. The function accepts raw count data from a temporal RNA-Seq dataset of two or more groups. The function outputs parameters mean, phase and amplitude are for each group.
#' @param countData	matrix containing non-negative integers; each column represents a sample, each row represents a gene/transcript.
#' @param group	vector containing the name of each sample.
#' @param time	vector containing numeric values of the time for each sample.
#' @param countData	matrix containing non-negative integers; each column represents a sample, each row represents a gene/transcript.
#' @param period	numeric value to indicate period length of the oscillation. Default for circadian data 24 h.
#' @param sample_name	vector containing sample names.
#' @param batch	vector containing potential batch effects between samples.
#' @param nthreads vector numeric value to indicate the threads for parallel computing .
#' @return a list that contains the following data.frames: results (summary of results), parameters (rhythmic parameters), ncounts (normalized counts), counts (raw counts), cook (cook's distance)
#' @examples countData = simData[["countData"]]
#' group = simData[["group"]]
#' time  = simData[["time"]]
#' dryList = dryseq(countData,group,time)
#' head(dryList[["results"]])    # data frame summarizing results
#' head(dryList[["parameters"]]) # coefficients: phase, amplitude and mean for each group
#' head(dryList[["ncounts"]])    # normalized counts
#' head(dryList[["counts"]])     # raw counts
#' head(dryList[["cook"]])       # cook's distance
#' @details DryR assesses rhythmicity and mean differences of gene expression in RNA-Seq count data.
#'     As proposed (Love et al. 2014), a count Y of a gene in a sample s can be modeled as a negative binomial with a fitted mean μ_gs and a gene-specific dispersion parameter θ_g.
#'     \cr \cr  \emph{Y_gs}~NB(\emph{μ_gs},\emph{θ_g})\cr\cr
#'     The fitted mean is proportional to the quantity q of fragments that correspond to a gene in a sample scaled by a sample-specific scaling factor s_s (Love et al. 2014).
#'     This scaling factor depends on the sampling depth of each library and can be estimated using the median-of-ratios method of DESeq2 (Anders and Huber 2010).
#'     \cr \cr   \emph{μ_gs} = \emph{s_s q_gs}\cr \cr
#'     DryR estimates gene-specific distribution θg using empirical Bayes shrinkage described by Love et al. (Love et al. 2014).
#'     and variance is computed from the following relationship to the dispersion parameter θ:
#'     \cr \cr   Var(\emph{Y_gs}) = E[(\emph{Y_gs} + \emph{θ_g} \emph{μ^2_gs})]\cr \cr
#'     The fit uses a generalized linear model with a logarithmic link function. Sample specific size factor (s_s) is defined as an offset. The full GLM is defined as follows:
#'     \cr \cr log2(\emph{μ_gct}) = \emph{m_g} + \emph{m_gc} + \emph{α_gc} cos(\emph{ωt}) + \emph{β_gc} sin(\emph{ωt}) + log2(\emph{s_s})\cr \cr
#'      μ is the raw count for gene g, condition/group c and Zeitgeber/circadian time t. α and β are coefficients of the cosine and sine functions, respectively. m is a coefficient to describe a mean expression level.
#'      When necessary, a batch specific mean (m) can be given to the dryseq function to account for technical batch effects.
#'      A technical batch effect is not allowed to be confounding so the resulting model matrix is fully ranked.
#'      To select an optimal gene-specific model, dryseq first assesses rhythmicity across the different conditions. To this end, dryR defines different models across all groups.
#'      Models refined to have either zero (non-rhythmic pattern) or non-zero (rhythmic pattern) α and β coefficients for each analyzed group. Moreover, for some models the values of α and β can be also shared within any combination of all groups
#'      The coefficients α and β were used to calculate the phase (arctan(α/β)) and amplitude (log2-fold change peak-to-trough; 2sqrt(α^2+β^2) ) of a gene.
#'      Bayesian information criterion (BIC) based model selection was employed to account for model complexity using the following formula:
#'      \cr \cr   BIC_j = ln(n)k - 2ln(L̂)\cr \cr
#'      L̂ is defined as the log-likelihood of the model j from the regression, n is the number of data points and k is the number of parameters.
#'      To assess the confidence of the selected model j we calculated the Schwarz weight (BICW):
#'      \cr \cr   BICW_j = e^(0.5ΔBIC_j)\ sum(e^0.5 ΔBIC_m), with ΔBIC_j - BIC_j - BIC_m*\cr \cr
#'      m* is the minimum BIC value in the entire model set. Dryseq consideres the BICW_j as the confidence level for model j. The model with the highest BICW is selected as the optimal model within the set of all defined models.
#'      In a second iteration step, dryR set the coefficient α and β to the values of the selected model in the first regression.
#'      dryseq then defined different models for the mean coefficient with differing or shared means between groups. Each model is solved using generalized linear regression and each gene was assigned to a preferred model based on the BICW as described above for the first iteration.
#'      The model selection is sensitive to outliers: dryseq provide a cook's distance for each gene. A fit for a gene with a maximum cook's distance of higher than 1 should be considered with care.
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology
#' @references Anders, S. and Huber, W. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology
dryseq=function(countData,group,time,period=24,sample_name=colnames(countData),batch=rep("A",length(sample_name)),n.cores=round(detectCores()*.6,0) ){
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

  s1 <- sin(2*pi*time/period)
  c1 <- cos(2*pi*time/period)

  conds  = cbind(group,s1,c1,batch)
  colnames(conds) = c("group","s1","c1","batch")

  colData <- data.frame(row.names=colnames(countData), conds)
  N=length(unique(group))

  ############################
  # FIT RHYTHMS

  message("fitting rhythmic models")

  models = create_matrix_list(time,N,period)
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

  dds = DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData,   design = models[[length(models)]])
  dds.full = DESeq2::DESeq(dds, full=models[[length(models)]], betaPrior = F, fitType = "parametric", test = "Wald", parallel =T, quiet = T)

  deviances = sapply(models[-length(models)], function(m){
    dds.x = DESeq2::nbinomWaldTest(dds.full, modelMatrix= m, betaPrior = F, quiet = T)
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
        dds.m <- DESeq2::nbinomWaldTest(dds.m[gene], modelMatrix= as.matrix(m), betaPrior = F) # Re-run wald test
        #return(list(dds.m,mcols(dds.m)$deviance)) # Returning deviances (-2 * log likelihood) // https://support.bioconductor.org/p/107472/
        return(list(dds.m, mcols(dds.m)$deviance)) # Returning deviances (-2 * log likelihood) // https://support.bioconductor.org/p/107472/

      })
    }

    if(length(gene)==0){dev = list (NA, NA)}



    return(dev)
  }

  deviance_mean = NULL
  for (cm_r in 1:length(models)){

    if(!is.na(DDS_dev[[cm_r]][1])){
      deviance_mean.x  = rbind(sapply(1:15,function(x) {DDS_dev[[cm_r]][[x]][[2]]}))
      rownames(deviance_mean.x)  = rownames(DDS_dev[[cm_r]][[1]][[1]])
      deviance_mean    = rbind(deviance_mean, deviance_mean.x)}
  }

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
    out = compute_param(dds, gene ,period,N)
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

  # to add flags for low expression (counts), high cook's distance
  # to be added error messages when only one group is given etc.
}
