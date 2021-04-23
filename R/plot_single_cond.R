#' This function plots data and fit for one condition
#' @export
plot_single_cond=function(out,gene_name){
  out.ncount = out$results[,1:(ncol(out$results)-7)]
  x=out.ncount[gene_name,]
  period=out$period
  t=out$time
  t.2=seq(min(t),max(t),0.1)
  c1=out$results[gene_name,'c1']
  s1=out$results[gene_name,'s1']

  u=out$results[gene_name,'Intercept']
  c1.f=c1*cos(2*pi*t.2/period)
  s1.f=s1*sin(2*pi*t.2/period)

  fit=u+c1.f+s1.f
  df=data.frame(t=t,dat=log2(1+as.numeric(x)))
  df.2=data.frame(t.2=t.2,fit=fit)
  g1=ggplot(data=df,aes(x=t,y=dat)) +
    geom_point(shape=21,color='white',fill='black',size=2) +
    geom_line(data=df.2,aes(x=t.2,y=fit),alpha=.8) +
    theme_bw() + theme(aspect.ratio=1,axis.text=element_text(size = 12),plot.title = element_text(size=10)) +
    xlab("Time") + ylab("Normalized read count log2(1+x)")  +
    ggtitle(paste(gene_name,"\n","adjpval:",format(out$results[gene_name,'padj'],scientific=T,digits = 3),", phase:",round(out$results[gene_name,'phase'],2),
                  ", amp:",round(out$results[gene_name,'amp'],2)))
  print(g1)
}
