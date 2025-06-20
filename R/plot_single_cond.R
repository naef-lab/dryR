#' This function plots data and fit for one condition for glms
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
    xlab("Time") + ylab("Log2 normalized counts") +
    ggtitle(paste(gene_name,"\n","adj_pval:",format(out$results[gene_name,'padj'],scientific=T,digits = 3),", phase:",round(out$results[gene_name,'phase'],2),
                  ", amp:",round(out$results[gene_name,'amp'],2)))
  print(g1)
}

#' This function plots data and fit for one condition for linear models
#' @export
plot_single_cond_lm <- function(out, gene_name) {
  # extract normalized counts (all columns except last 9 stats)
  counts <- out$results[, seq_len(ncol(out$results) - 9)]
  expr   <- counts[gene_name, ]
  period <- out$period
  times  <- out$time

  # fine grid for smooth fit
  times_grid <- seq(min(times), max(times), by = 0.1)

  # coefficients from cycler
  intercepts <- out$results[gene_name, grep("Intercept|mean", colnames(out$results))]
  coef_cos   <- out$results[gene_name, "c1"]
  coef_sin   <- out$results[gene_name, "s1"]

  # fitted cos/sin curve
  cos_part <- coef_cos * cos(2 * pi * times_grid / period)
  sin_part <- coef_sin * sin(2 * pi * times_grid / period)
  fit_vals <- intercepts + cos_part + sin_part

  # data frames for plotting
  df_points <- data.frame(time = times, expression = as.numeric(expr))
  df_curve  <- data.frame(time = times_grid, fit = fit_vals)

  # build plot
  p <- ggplot(df_points, aes(x = time, y = expression)) +
    geom_point(shape = 21, color = "white", fill = "black", size = 2) +
    geom_line(data = df_curve, aes(x = time, y = fit), alpha = 0.8) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      axis.text    = element_text(size = 12),
      plot.title   = element_text(size = 10)
    ) +
    labs(
      x = "Time",
      y = "Log2 normalized counts",
      title = sprintf(
        "%s\nadj_pval: %.3g, phase: %.2f, amp: %.2f",
        gene_name,
        out$results[gene_name, "padj"],
        round(out$results[gene_name, "phase"], 2),
        round(out$results[gene_name, "amp"], 2)
      )
    )

  print(p)
}


