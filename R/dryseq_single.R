#' @export
dryseq_single=function(countData, 
                sample_name=names(countData),
                group, 
                single,
                time,
                period=24){

  sel =  group %in% single
  
  if(!any(sel)){warning("Your single condition is not included in 'group'")}
  
  time      = time[sel]
  group     = group[sel]
  countData = countData[,sel]
  sample_name = sample_name[sel]
  
  countData = countData[rowSums(countData)!=0,]
  
  s1 <- sin(2*pi*time/period)
  c1 <- cos(2*pi*time/period)
  
  conds  = cbind(s1,c1)
  colnames(conds) = c("s1","c1")
  
  colData <- data.frame(row.names=colnames(countData), conds)
  N=length(unique(group))
  
  ############################
  # FIT RHYTHMS
  ######################
  
  dds = DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ s1 + c1)  
  dds = DESeq(dds, test="LRT", reduced=~1)
  
  ################
  # coefficients / mean, amplitude, phase, p-value
  ###############
  
  res = as.data.frame(cbind(results(dds),coefficients(dds)))
  res = res[,c('pvalue','padj','Intercept','s1','c1')]
  
  phase=period/(2*pi)*atan2(res$s1,res$c1)
  phase=phase%%period
  amp =2*sqrt(res$s1^2+res$s1^2)
  
  res=data.frame(res,phase,amp)
  
  
  #normalized counts
  ncounts       = counts(dds, normalized = TRUE)
  global_df = as.data.frame(cbind(ncounts,res))

  
  out = list()
  
  out[["time"]]  = time
  out[["period"]] = period
  out[["results"]] = global_df
  out[["single"]] = single
  message("finished!")
  return(out)
  
}
