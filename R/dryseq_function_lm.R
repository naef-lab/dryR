#' Differential rhythmicity analysis for RNA-Seq datasets
#'
#' This function performs a rhythmicity analysis based on linear models with a subsequent models selection. The function accepts raw count data from a temporal RNA-Seq dataset of two or more groups. The function outputs parameters mean, phase and amplitude are for each group.
#' @param countData	matrix containing library normalized and log transformed reads counts; each column represents a sample, each row represents a gene/transcript.
#' @param group	vector containing the name of each sample.
#' @param time	vector containing numeric values of the time for each sample.
#' @param period	numeric value to indicate period length of the oscillation. Default: circadian data period of 24 h.
#' @param sample_name	vector containing sample names. Default: colnames are sample names.
#' @param batch	vector containing potential batch effects between samples. Default: no batch effect.
#' @param nthreads vector numeric value to indicate the threads for parallel computing .Default: 60 \% of detected cores.
#' @return a list that contains the following data.frames: results (summary of results), parameters (rhythmic parameters), ncounts (normalized counts), counts (raw counts), cook (cook's distance)
#' @examples countData = simData[["countData"]]
#' group = simData[["group"]]
#' time  = simData[["time"]]
#' dryList = dryseq(countData,group,time)
#' head(dryList[["results"]])    # data frame summarizing results
#' head(dryList[["parameters"]]) # coefficients: phase, amplitude and mean for each group
#' head(dryList[["ncounts"]])    # normalized counts
#' @details DryR assesses rhythmicity and mean differences of gene expression in RNA-Seq count data.
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

dryseq_lm=function(countData,group,time,period=24,sample_name=colnames(countData),batch=rep("A",length(sample_name)),n.cores=round(detectCores()*.6,0) ){
  require("DESeq2")
  require("combinat")
  require("parallel")
  require("gplots")
  library("doMC")
  library("circular")
  library("RColorBrewer")

  sel       = order(group,time)
  time      = time[sel]
  group     = group[sel]
  countData = countData[,sel]


  register(MulticoreParam(n.cores))
  registerDoMC(cores=n.cores)

  #countData = countData[rowSums(countData)!=0,]

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
  models = lapply(models, function(l) l[,c(grep("u",colnames(l)),grep("a|b",colnames(l)))])

  for (i in 1:length(models)){
    rownames(models[[i]]) = rownames(colData)}

  if (length(unique(batch))>1) {
    # add the batch effect
    model_b = as.matrix(model.matrix(~  batch),contrasts.arg=NULL)[,2:length(unique(batch)),drop=F]
    colnames(model_b)=paste0("BATCH_",unique(batch)[-1])
    models = lapply(models, function(l) cbind(model_b,l))
    models = lapply(models, function(l) l[,c(grep("u",colnames(l)),grep("BATCH",colnames(l)),grep("a|b",colnames(l)))]   )
  }

  fit = parallel::mclapply(split(countData, rownames(countData)),
                           do_all_lm,
                           my_mat= models,
                           mc.cores=n.cores)

  # calculate the BIC
  BIC = unlist(fit)[grep('BIC$',names(unlist(fit)))]
  BIC= matrix(BIC,nrow=nrow(countData),byrow=T)

  #calculate the BICW
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

  gene.list=as.list(rownames(countData))
  fit = parallel::mclapply(gene.list,
                           FUN=do_all_lm_mr,
                           countData,
                           my_mat_r = models,
                           my_mat_m = model_mean_cond,
                           choosen_model = choosen_model,
                           mc.cores=n.cores)

  #extract BIC
  BIC_mean = unlist(fit)[grep('BIC$',names(unlist(fit)))]
  BIC_mean = matrix(BIC_mean,nrow=nrow(countData),byrow=T)

  #calculate the BICW
  BICW_mean = t(apply(BIC_mean,1,compute_BICW))

  choosen_model_mean = apply(BIC_mean,1,which.min)
  choosen_model_mean_BICW = apply(BICW_mean,1,max)

  ################
  # coefficients / mean, amplitude and phase
  ###############

  message("extracting rhythmic parameters")
  parameters=NULL

  parameters =  foreach (i = 1:nrow(countData)) %dopar% {
    gene = rownames(countData)[i]

    dds= fit[[i]][[choosen_model_mean[i]]]$param
    out = compute_param_l(dds,period, N)
    return(out)
  }
  parameters            = data.frame(t(do.call(cbind, parameters)))
  colnames(parameters)  = c(paste(c('mean','a','b','amp','relamp','phase'),rep(unique(group),each =6), sep = "_"))
  rownames(parameters)  = rownames(countData)


  #normalized counts
  ncounts_RF = countData

  # generate a table summarizing the analysis
  complete_parameters = cbind(parameters,choosen_model,choosen_model_BICW, choosen_model_mean, choosen_model_mean_BICW)
  global_table_df = merge(ncounts_RF,complete_parameters, by="row.names")

  rownames(global_table_df) = global_table_df[,1]
  global_table_df           = global_table_df[,-1]

  out = list()

  out[["results"]]     = global_table_df
  out[["BICW_rhythm"]] = BICW
  out[["BICW_mean"]]   = BICW_mean
  out[["values"]]      = ncounts_RF
  out[["parameters"]]  = complete_parameters

  message("finished!")
  return(out)


}
