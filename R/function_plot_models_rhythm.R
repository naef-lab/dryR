#' Visualization of dryseq results
#'
#' This function allows to plot the results of dryseq.
#' @param dryList from dryseq output
#' @param file_path_name folder to store output
#' @export pdf that summarize the results of dryseq
#' @examples
#' XYZn
plot_models_rhythm = function(dryList,file_path_name, period=24){
  if (length(dev.list()!=0)) {dev.off()}
  t=dryList$time
  group=dryList$group
  x=dryList$results
  
  if("vsd" %in% names(dryList)){
    x[,1:ncol(dryList$vsd)]=dryList$vsd
    
  }
  
  pdf(file =paste0(file_path_name,'summary_heatmap_models.pdf'), paper = "a4")
  nb = table(x[,'chosen_model'])
  nb = nb[order(-nb)]
  mo = as.numeric(names(nb))
  
  for(i in mo){
    
    x_s = subset(x, chosen_model==i)
    
    if(nrow(x_s)>1){
      
      pos_phase = grep('phase',names(x_s))
      
      sum_phase = apply(x_s[,pos_phase],2,sum,na.rm=T)
      if(sum(sum_phase)!=0){
        x_s = x_s[order(x_s[,pos_phase[min(which(sum_phase !=0))]]),]
      }
      
      x_s=as.matrix(x_s[,1:length(t)])
      
      ## Average replicates
      com= paste(t,group)
      com.u=unique(com)
      t.u=sapply(strsplit(com.u," "),"[[",1)
      cond.u=sapply(strsplit(com.u," "),"[[",2)
      
      ss=split(1:length(com),com)
      
      x_s.m=sapply(ss,function(x) if(length(x) > 1){rowMeans(x_s[,x])}else{x_s[,x]})
      
      #remove mean per condition
      condi=sapply(strsplit(colnames(x_s.m)," "),"[[",2)
      ss=split(1:length(condi),condi)
      
      for(k in ss){ x_s.m[,k]= sweep(x_s.m[,k],1,rowMeans(x_s.m[,k]),FUN="-") }
      x_s.m=x_s.m[,match(com.u,colnames(x_s.m))]
      
      
      pheatmap(as.matrix(x_s.m),
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               show_rownames = F ,
               scale = "row",
               gaps_col=which(diff(as.numeric(as.factor(cond.u)))!=0),
               main = paste("model",i," #Genes",nrow(x_s),sep =" "),
               col = colorRampPalette(c('blue','yellow'))(1000))
      
      x_s = subset(x, chosen_model==i)
      ba = 1
      gg = list()
      
      for(kk in pos_phase){
        
        if(sum(x_s[,kk],na.rm=T) !=0){
          gg[[ba]] = circular_phase24H_histogram(x_s[,kk], unique(condi)[ba], period)
          ba = ba + 1
        }
      }
      
      if(length(gg)>0){gridExtra::grid.arrange(grobs = gg, ncol = 3, nrow = 1+round(length(pos_phase)/3))}
      
      bas = unique(sum_phase)
      bas=bas[bas!=0]
      la = match(bas, sum_phase)
      if (length(la) > 1) {
        pairs(x_s[, pos_phase[la]], cex = 0.5)
      }
      
    }
    
  }
  dev.off()
  
}
