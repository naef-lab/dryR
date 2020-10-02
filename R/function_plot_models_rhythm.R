#' Visualization of dryseq results
#'
#' This function allows to plot the results of dryseq.
#' @param dryList from dryseq output
#' @param file_path_name folder to store output
#' @export pdf that summarize the results of dryseq
#' @examples
#' XYZn
plot_models_rhythm = function(dryList,file_path_name){

  t=dryList$time
  group=dryList$group
  x=dryList$results

  pdf(file =paste0(file_path_name,'summary_heatmap_models.pdf'))

  nb = table(x[,'choosen_model'])
  nb = nb[order(-nb)]
  mo = as.numeric(names(nb))

  for(i in mo){

    x_s = subset(x, choosen_model==i)

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

      x_s = subset(x, choosen_model==i)
      par(mfrow=c(1+round(length(pos_phase)/3),3))
      ba = 1
      for(kk in pos_phase){

        if(sum(x_s[,kk],na.rm=T) !=0){
          circular_phase24H_histogram(x_s[,kk], unique(condi)[ba], period)
        }
        ba = ba + 1
      }
      bas = unique(sum_phase)
      la = match(bas,sum_phase)
      if(length(unique(sum_phase[sum_phase!=0])) > 1){

        pairs(x_s[,pos_phase[la]],cex=0.5)

      }

    }

  }
  dev.off()

}


