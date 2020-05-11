#' Visualization of dryseq results
#'
#' This function allows to plot the results of dryseq.
#' @param DF results from dryseq
#' @param file_path_name folder to store output
#' @export pdf that summarize the results of dryseq
#' @examples
#' XYZ
plot_models_rhythm = function(DF,file_path_name,time,group, period=24){

  x=DF

  pdf(file =paste(file_path_name,'.pdf',sep=''))

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

      x_s=as.matrix(x_s[,1:length(time)])

      ## Average replicates
      com= paste(time,group)
      com.u=unique(com)
      t.u=sapply(strsplit(com.u," "),"[[",1)
      cond.u=sapply(strsplit(com.u," "),"[[",2)

      ss=split(1:length(com),com)

      x_s.m=sapply(ss,function(x) rowMeans(x_s[,x]))

      #remove mean per condition
      condi=sapply(strsplit(colnames(x_s.m)," "),"[[",2)
      ss=split(1:length(condi),condi)

      for(k in ss){ x_s.m[,k]= sweep(x_s.m[,k],1,rowMeans(x_s.m[,k]),FUN="/") }
      x_s.m=x_s.m[,match(com.u,colnames(x_s.m))]


      heatmap.2(as.matrix(x_s.m),
                Rowv = NA,
                Colv = NA,
                ylab = NA ,
                dendrogram = "none",
                labCol = paste('ZT',t.u, sep = "_"),
                labRow = NA ,
                scale = NULL,
                trace='none',
                colsep=which(diff(as.numeric(as.factor(cond.u)))!=0),
                main = paste("model",i," #Genes",nrow(x_s),sep =" "),
                col = colorRampPalette(c('blue','black','yellow'))(1000),
                key=F)


      par(mfrow=c(1+round(length(pos_phase)/3),3))
      ba = 1
      for(kk in pos_phase){

        if(sum(x[,kk],na.rm=T) !=0){
          circular_phase24H_histogram(x[,kk], unique(condi)[ba], period)
        }
        ba = ba + 1
      }
      bas = unique(sum_phase)
      la = match(bas,sum_phase)
      if(length(unique(sum_phase)) > 1){

        pairs(x[,pos_phase[la]],cex=0.5)

      }

    }

  }
  dev.off()

}


