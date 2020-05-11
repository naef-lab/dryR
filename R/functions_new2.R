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

simply_it.2 = function(x){

  a = match(x,x)
}


#####################################
create_matrix_list = function(t,co,period){
  my_matrix = list()

  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)

  MAT = cbind(rep(1,length(t)/co),c[1:(length(t)/co)],s[1:(length(t)/co)])
  GMAT = matrix(NA,ncol=3*co, nrow =length(t))
  colnames(GMAT) = c(paste(c('u','a','b'),rep(1:co,each =3), sep = "."))
  for(i in 1:co){
    border = length(t)/co
    s = (i-1)*border + 1
    e = i*border
    sr = (i-1)*3 + 1
    er = i*3
    GMAT[s:e,sr:er] = MAT
  }


  sum_m = 0
  for(i in 0:co){
    sum_m = sum_m + comb(co,i)
  }
  vn = rep(F,co)
  for(i in 1:co){
    g = rep(F,co)
    g[1:i] = TRUE

    p = unique(permn(g))
    v = matrix(unlist(p),ncol = co,byrow = TRUE)
    vn = rbind(vn,v)

  }


  mp = seq(1,3*co,3)
  vn = vn[,rep(1:co,each=3)]
  vn[,seq(1,3*co,3)] = TRUE
  vn = data.frame(vn)
  vn[,dim(vn)[2] + 1]=(apply(vn,1,nbt)-co)/2
  colnames(vn) = c(paste(c('u','a','b'),rep(1:co,each =3), sep = "."),'nb_cycl')

  model = 1
  for(g in 0:co){


    nb_cycl =g
    com = expand.grid(rep(list(1:nb_cycl),nb_cycl))
    simply = apply(com,1,simply_it)
    poss =match(unique(simply),simply)
    com_l = com[poss,]
    pos = which(vn$nb_cycl == g)

    for(k in pos){
      if(g > 1){
        for(v in 1:nrow(com_l)){
          gmat = GMAT[,unlist(vn[k,-ncol(vn) ])]
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
              for( z in 1:((length(he)-2)/2)){
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
        gmat = GMAT[,unlist(vn[k,-ncol(vn) ])]
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

create_matrix_list_mean=function(N,group){
  com = expand.grid(rep(list(1:N),N))
  simply = as.data.frame(t(apply(com,1,simply_it.2)))
  simply = do.call("paste",simply)
  poss =match(unique(simply),simply)
  com_l = com[poss,]
  names(com_l)=unique(group)
  com_l=com_l[order(apply(com_l,1,function(x) length(unique(x))),apply(com_l,1,function(x) length(which(x==max(x))))),]
  rownames(com_l)=1:nrow(com_l)
  heatmap.2(as.matrix(com_l),
            dendrogram='none',
            Rowv=FALSE,
            Colv=FALSE,
            trace='none',
            col=brewer.pal(n = 8, name = 'Set2'),
            colsep=1:ncol(com_l),
            rowsep=1:nrow(com_l),
            cexRow=0.8,
            cexCol = 0.8,
            key = FALSE)


  com_l=com_l[,rep(1:N,times=table(group))]
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


annotate_matrix=function(m,group){
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

############################


compute_BICW=function(x){
  x = as.numeric(x)
  BIC_min = min(x)
  test = exp(-0.5*(x-BIC_min))/sum(exp(-0.5*(x-BIC_min)))
  return(test)
}

#########################

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



########################

dryseq=function(countData,group,time,T_=24,sample_name=colnames(countData),batch=rep("A",length(sample_name)),n.cores=round(detectCores()*.6,0) ){

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

#####################################
make_circ_coord = function(t,x,ttot) {
  dt=(t[2]-t[1])*.45
  a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
  h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
  list(angles=a,heights=h)
}
#####################################
circular_phase24H_histogram<-function(x,name,ttot){
  color_hist = rgb(0.6,0,0.2)
  br=0:ttot
  h=hist(x, br=br,plot=F)
  co=make_circ_coord(br[-1],h$counts,ttot)
  radial.plot(co$heights,co$angle,br[-1]-br[2],
              clockwise=T,start=pi/2,main=paste("",name),
              rp.type='p',poly.col=color_hist, xlab = "",ylab = "", show.grid.labels=0)
}
#################################

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


