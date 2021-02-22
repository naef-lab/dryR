#####################################
comb = function(n,k){
  factorial(n)/(factorial(k)*factorial(n-k))
}
#####################################
nbt = function(x){
  l=length(which(x) == TRUE)
  l
}
#####################################
simply_it = function(x){
  a = 0
  for(i in x) {
    a= paste(a,paste(which(x == as.numeric(i)), collapse = "",sep = ""), collapse = "", sep = "")
  }
  a
}
#####################################
simply_it.2 = function(x){

  a = match(x,x)
}
#####################################
make_circ_coord = function(t,x,ttot) {
  dt=(t[2]-t[1])*.45
  a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
  h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
  list(angles=a,heights=h)
}
#####################################
circular_phase24H_histogram = function(x,name,ttot){
  color_hist = rgb(0.6,0,0.2)
  br=0:ttot
  h=hist(x, br=br,plot=F)
  co=make_circ_coord(br[-1],h$counts,ttot)
  radial.plot(co$heights,co$angle,br[-1]-br[2],
              clockwise=T,start=pi/2,main=paste("",name),
              rp.type='p',poly.col=color_hist, xlab = "",ylab = "", show.grid.labels=0)
}
#####################################
compute_BICW = function(x){
  x = as.numeric(x)
  BIC_min = min(x)
  test = exp(-0.5*(x-BIC_min))/sum(exp(-0.5*(x-BIC_min)))
  return(test)
}
#################################3
compute_RSS = function(x, matX){

  xx = solve(t(matX)%*%matX)
  y = xx %*% t(matX) %*% as.numeric(x)
  y = as.matrix(y)
  rownames(y)  = colnames(matX)
  RSS = t(x) %*% x -t(x) %*% matX %*% xx %*% t(matX) %*% x
  list(param=y,RSS = RSS)
}

#############################################
compute_BIC = function(A,n){

  p = length(A$param)
  #AIC = n * log(A$RSS/n, base = exp(1)) + 2* p + 2*p*(p +1) /(n-p-1)
  BIC=  n * log(A$RSS/n, base = exp(1))  + log(n, base = exp(1)) * p
  list(BIC = BIC, param = A$param)

}

##############################
do_all_lm = function(x,my_mat){
  x = as.numeric(x)
  n = length(x)

  my_fit = lapply(my_mat,compute_RSS, x = x)
  my_BIC =lapply(my_fit,compute_BIC,n=n)
  my_BIC
}

##########################

do_all_lm_mr = function(x,countData,my_mat_r, my_mat_m, choosen_model){
  i=match(x,rownames(countData))
  x=countData[x,]
  M=my_mat_r[[choosen_model[i]]]
  #build the gene specific model from the rhythmic point of view
  gene_specific_mean_models = lapply(my_mat_m,
                                     function(x) cbind(x,M[,-grep("u",colnames(M))]))
  x = as.numeric(x)
  n = length(x)

  my_fit = lapply(gene_specific_mean_models,compute_RSS, x = x)
  my_BIC =lapply(my_fit,compute_BIC,n=n)
  my_BIC
}

####################################
compute_param = function(dds, gene, period=T_,N){

  dds = dds[gene,]
  param = c(paste(rep(c('u','a','b'),each=N),rep(1:N,3), sep = "."))

  paramout = rep(NA,N*6)

  for(i in 1:N){

    u=coef(dds)[grep(paste(param[i],"Intercept",sep="|"), colnames(coef(dds)))]
    a=coef(dds)[grep(param[i+N], colnames(coef(dds)))]
    b=coef(dds)[grep(param[i+N*2], colnames(coef(dds)))]

    if(length(u) ==0) u=NA
    if(length(a) ==0) a=NA
    if(length(b) ==0) b=NA

    phase=period/(2*pi)*atan2(b,a)
    amp =2*sqrt(a^2+b^2)
    relamp=0.5*amp/u
    if(!is.na(phase)){
      #if(phase<0) phase=phase+period
      #if(phase>period) phase=phase-period
      phase=phase%%period
    }
    paramout[(1:6 + 6*(i-1))] = c(u,a,b,amp,relamp,phase)
  }

  #names(paramout) = c(paste(c('mean','a','b','amp','relamp','phase'),rep(1:N,each =6), sep = "_"))
  paramout
}
####################################
compute_param_l = function(dds, period=T_, N){


  param = c(paste(rep(c('u','a','b'),each=N),rep(1:N,3), sep = "."))

  paramout = rep(NA,N*6)

  for(i in 1:N){

    u=dds[grep(paste(param[i],"Intercept",sep="|"), rownames(dds)),1]
    a=dds[grep(param[i+N], rownames(dds)),1]
    b=dds[grep(param[i+N*2], rownames(dds)),1]

    if(length(u) ==0) u=NA
    if(length(a) ==0) a=NA
    if(length(b) ==0) b=NA

    phase=period/(2*pi)*atan2(b,a)
    amp =2*sqrt(a^2+b^2)
    relamp=0.5*amp/u
    if(!is.na(phase)){
      #if(phase<0) phase=phase+period
      #if(phase>period) phase=phase-period
      phase=phase%%period
    }
    paramout[(1:6 + 6*(i-1))] = c(u,a,b,amp,relamp,phase)
  }

  #names(paramout) = c(paste(c('mean','a','b','amp','relamp','phase'),rep(1:N,each =6), sep = "_"))
  paramout
}

#####################################
create_matrix_list = function(t, conds, n.co, period){
  require(combinat)

  my_matrix = list()

  c <- cos(2*pi*t/period)
  s <- sin(2*pi*t/period)

  MAT <- cbind(rep(1,length(t)),c[1:length(t)],s[1:length(t)])
  GMAT <- matrix(NA,ncol=3*n.co, nrow =length(t))
  rownames(GMAT) <- conds
  colnames(GMAT) <- c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."))

  it <- 1
  for(i in unique(rownames(GMAT))){
    GMAT[rownames(GMAT)==i,grep(paste0('.',it,'$'),colnames(GMAT))] = MAT[rownames(GMAT)==i,]
    it=it+1
  }

  vn = rep(F,n.co)
  for(i in 1:n.co){
    g = rep(F,n.co)
    g[1:i] = TRUE
    p = unique(combinat::permn(g))
    v = matrix(unlist(p),ncol = n.co,byrow = TRUE)
    vn = rbind(vn,v)

  }


  vn = vn[,rep(1:n.co,each=3)]
  vn[,seq(1,3*n.co,3)] = TRUE
  vn = data.frame(vn,row.names= NULL)
  vn[,dim(vn)[2] + 1]=(apply(vn,1,nbt)-n.co)/2
  colnames(vn) = c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."),'nb_cycl')

  model = 1
  for(g in 0:n.co){


    nb_cycl =g
    com = expand.grid(rep(list(1:nb_cycl),nb_cycl))
    simply = apply(com,1,simply_it)
    poss =match(unique(simply),simply)
    com_l = com[poss,]
    pos = which(vn$nb_cycl == g)

    for(k in pos){
      if(g > 1){
        for(v in 1:nrow(com_l)){
          gmat = GMAT[,unlist(vn[k,-dim(vn)[2]])]
          ve = as.numeric(com_l[v,])
          id =1
          sa = ve
          while(length(ve) !=0){

            poc = which(sa == ve[1])
            po = which(ve ==ve[1])
            if(length(poc) !=1){
              poch =c(2*poc-1,2*poc)
              poch =poch[order(poch)]
              he = grep("[ab]",colnames(gmat))
              he = he[poch]
              pp=0
              for(z in 1:((length(he)-2)/2)){
                repl1 = which(gmat[,he[2*z+1]]!='NA')
                repl2 = which(gmat[,he[2*z+2]]!='NA')
                gmat[repl1,he[1]] = gmat[repl1,he[2*z+1]]
                gmat[repl2,he[2]] = gmat[repl2,he[2*z+2]]
                colnames(gmat)[he[1]]= paste(colnames(gmat)[he[1]],colnames(gmat)[he[2*z+1]],sep=',')
                colnames(gmat)[he[2]]= paste(colnames(gmat)[he[2]],colnames(gmat)[he[2*z+2]],sep=',')
                gmat[repl1,he[2*z+1]] =NA
                gmat[repl2,he[2*z+2]]=NA
                pp = pp+2
              }
              id = id+1
              ve = ve[-po]
            }else{
              ve = ve[-1]
            }

          }
          gmat[is.na(gmat)] =0
          del=which(apply(gmat,2,function(x) length(which(x == 0))) == length(t))
          if(length(del)!=0){
            gmat = gmat[,-del]
          }
          my_matrix[[model]] = gmat
          model = model + 1
        }
      }else{
        gmat = GMAT[,unlist(vn[k,-dim(vn)[2]])]
        gmat[is.na(gmat)] =0
        del =which(apply(gmat,2,function(x) length(which(x == 0))) == length(t))
        if(length(del)!=0){
          gmat = gmat[,-del]
        }
        my_matrix[[model]] = gmat
        model = model +1
      }
    }

  }



  return(my_matrix)
}
#####################################
create_matrix_list_mean = function(N,group){
  com = expand.grid(rep(list(1:N),N))
  simply = as.data.frame(t(apply(com,1,simply_it.2)))
  simply = do.call("paste",simply)
  poss =match(unique(simply),simply)
  com_l = com[poss,]
  names(com_l)=unique(group)
  com_l=com_l[order(apply(com_l,1,function(x) length(unique(x))),apply(com_l,1,function(x) length(which(x==max(x))))),]
  rownames(com_l)=1:nrow(com_l)

  com_l=com_l[,match(group,names(com_l))]
  p=list()
  for(j in 1:nrow(com_l)){
    if(j==1){
      p[[j]] = as.matrix(rep(1,ncol(com_l)))

    }else{
      p[[j]]= model.matrix(~0+ factor(as.numeric(com_l[j,])))
    }
  }
  p
}
#####################################
annotate_matrix = function(m,group){
  if(ncol(m)==1){
    colnames(m)=paste("u",1:length(unique(group)),sep=".",collapse=".")
  }else{
    pos_ind= match(unique(group),group)
    m=as.matrix(m)
    l=list()
    for(k in 1:ncol(m)){
      l[[k]]=as.numeric(which(m[pos_ind,k]==1))
    }
    colnames(m)=sapply(l,function(x) paste("u",x,sep=".",collapse="."))
  }
  m
}

#####################################
dry_plot = function (dryList, gene)
{


  normal = FALSE
  if("ncounts" %in% names(dryList)){vsd        = log2(dryList[["ncounts"]]+1)}
  if("values" %in% names(dryList)){vsd         = dryList[["values"]]
                                    normal = T}

  parameters = dryList[["parameters"]][,grep("^mean|^a_|^b_|^amp|^phase|^relamp",colnames(dryList[["parameters"]]))]

  ID = rownames(dryList[["results"]] )[grep(paste0('^',gene,'$'),rownames(dryList[["results"]] ))]

  print(ID)

  d = vsd[ID, ]
  d = reshape2::melt(d)

  d$group            = dryList[["group"]]

  d$time            = as.numeric(dryList[["time"]])
  d$time            = d$time%%24

  suppressWarnings({ d <- Rmisc::summarySE(d, measurevar="value", groupvars=c("time","group")) })

  v = seq(0,24,1)
  fit_d_0 = parameters[which(rownames(parameters)==ID),grep("mean",colnames(parameters))] # intercept
  fit_d_1 = parameters[which(rownames(parameters)==ID),grep("a_",colnames(parameters))] # coefficient a
  fit_d_2 = parameters[which(rownames(parameters)==ID),grep("^b_",colnames(parameters))] # coefficient b

  fit_d_0[is.na(fit_d_0)] = 0
  fit_d_1[is.na(fit_d_1)] = 0
  fit_d_2[is.na(fit_d_2)] = 0

  m = data.frame(v)

  dd = data.frame(v)
  dd$v = v

  fit_values = function (x,n)
  { as.numeric((fit_d_0[n] + fit_d_1[n]*cos(2*pi*x/24)  + fit_d_2[n]*sin(2*pi*x/24)))  }


  for (u in 1:length(unique(d$group))){
    m[,u+1]  = NA
    m[,u+1]  = apply(dd,1, fit_values,u)
  }

  m = m[,-1]

  colnames(m) =  unique(dryList[["group"]])

  m =  reshape2::melt(m)
  m$time = rep(v, length(unique(d$group)))

  colnames(m)       = c("group","value","time")

  if(normal==FALSE) {m$value[which(m$value<0)] = 0}

  gg1 = ggplot2::ggplot(d, aes(x=time, y=value, group=group, color=group)) +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.4) +
    geom_point(size=2, shape=19) +
    xlab("Time (h)") +
    ggtitle(ID) +
    scale_x_continuous(breaks=c(0,6,12,18,24,30)) +
    theme_bw(base_size = 10) +
    theme(aspect.ratio = 1, panel.grid.minor=element_blank(), legend.position = "right") +
    geom_line(aes(x=time, y=(value), group=group), data = m, position=position_dodge(width=0.5)) +
    facet_wrap(~group)

  print(gg1)
}

