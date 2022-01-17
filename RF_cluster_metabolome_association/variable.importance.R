subset.bprf=function(bprfObj, idx.tree){
  subset=list()
  obj=bprfObj
  ntree=length(obj$forest)
  for(i in 1:length(obj)){
    if(length(obj[[i]])!=ntree) subset[[i]]=obj[[i]] else subset[[i]]=obj[[i]][idx.tree]
  }
  names(subset)=names(obj)
  return(subset)
}

importance.one=function(bprfObj, permutation=100, idx.x, nThreads=1){
  ### Find subset trees
  idx.tree=which(sapply(bprfObj$xinB, function(x)as.numeric(idx.x %in% x))==1)  
  x.obj=subset(bprfObj, idx.tree=idx.tree)
  Xtest=x.obj$Xtest
  forest=x.obj$forest
  xinB=x.obj$xinB
  ooB=x.obj$ooB
  Yhat=x.obj$Yhat[,idx.tree,drop=F]
  Yhat.fin=apply(Yhat, 1, function(yy)mean(yy[yy!=0]))
  Ytest.fin=x.obj$Y
  X.block=x.obj$X.block
  x.names=colnames(X.block)[idx.x]
  nn=length(Ytest.fin)
  package=x.obj$package
  
  if(sum(is.na(Yhat.fin))>0){
    Ytest.fin=x.obj$Ytest[!is.na(Yhat.fin),,drop=F]
    Yhat.fin=matrix(Yhat.fin[!is.na(Yhat.fin)])
  }
  
  r2.forest0=summary(lm(Ytest.fin~Yhat.fin))$r.squared
  rmse.forest0=sqrt(sum((Ytest.fin-Yhat.fin)^2)/length(Ytest.fin))
  accuracy.forest0=mean(abs(Ytest.fin-Yhat.fin)/abs(Ytest.fin))
  
  r2.forest1=rmse.forest1=accuracy.forest1=rep(0, permutation)
  
  perm=function(X.block, x.names, ooB, xinB, Ytest.fin, forest){
    idx.pm=sample(1:nrow(X.block),nrow(X.block),replace=F)
    xb1=X.block; xb1[,x.names]=xb1[idx.pm, x.names]
    Xtest1=lapply(1:length(ooB), function(i)xb1[ooB[[i]],xinB[[i]],drop=F])
    Yhat1=matrix(0, length(Ytest.fin), length(forest))
    for(i in 1:length(forest))
      Yhat1[ooB[[i]],i]= if(package=="ranger") predict(forest[[i]], data=data.frame(Xtest1[[i]]))$predictions else predict(forest[[i]], newdata=Xtest1[[i]])
    
    Yhat.fin1=apply(Yhat1, 1, function(y)mean(y[y!=0]))
    r2=summary(lm(Ytest.fin~Yhat.fin1))$r.squared
    rmse=sqrt(sum((Ytest.fin-Yhat.fin1)^2)/length(Ytest.fin))
    accuracy=mean(abs(Ytest.fin-Yhat.fin1)/abs(Ytest.fin))
    
    return(c(r2=r2, rmse=rmse, accuracy=accuracy))  
}

  if(nThreads==1){
    pb=txtProgressBar(min=1, max=permutation, style=3, char="*")
    for(p in 1:permutation){
      pm.data=perm(X.block=X.block, x.names=x.names, ooB=ooB, xinB=xinB,Ytest.fin=Ytest.fin, forest=forest)
      r2.forest1[p]=pm.data["r2"]
      rmse.forest1[p]=pm.data["rmse"]
      accuracy.forest1[p]=pm.data["accuracy"]
      setTxtProgressBar(pb, p)
    }
  } else{
    cl=parallel::makeCluster(nThreads);doParallel::registerDoParallel(cl)
    eval=foreach(p=1:permutation, .combine=rbind, .packages="randomForest") %dopar% {
      perm(X.block=X.block, x.names=x.names, ooB=ooB, xinB=xinB,Ytest.fin=Ytest.fin, forest=forest)
    }
    stopImplicitCluster(); parallel::stopCluster(cl)
    r2.forest1=eval[,"r2"]
    rmse.forest1=eval[,"rmse"]
    accuracy.forest1=eval[,"accuracy"]
  }
    
  dR2=mean(r2.forest0-r2.forest1); dR2.pvalue=(sum(r2.forest1>=r2.forest0)+1)/(1+permutation)
  dRMSE=mean(rmse.forest1-rmse.forest0); dRMSE.pvalue=(sum(rmse.forest1<=rmse.forest0)+1)/(1+permutation)  
  dAccuracy=mean(accuracy.forest1-accuracy.forest0); dAccuracy.pvalue=(sum(accuracy.forest1<=accuracy.forest0)+1)/(1+permutation) 

  return(c(dR2=dR2, dR2.pvalue=dR2.pvalue, dRMSE=dRMSE, dRMSE.pvalue=dRMSE.pvalue,dAccuracy=dAccuracy, dAccuracy.pvalue=dAccuracy.pvalue))
}

  
importance=function(bprfObj, permutation=100, idx.x, verbose=1, padjust=T, nThreads=1){
  X.block=bprfObj$X.block
  x.names=colnames(X.block)[idx.x]
  res=NULL

  for(i in 1:length(idx.x)){
   
  res=rbind(res, importance.one(bprfObj, permutation=permutation, idx.x=idx.x[i], nThreads=nThreads))
  unit=ifelse(length(idx.x)>10,round(length(idx.x)/10),1)
    if(verbose>0 & i%%unit==0){
        if(i==unit) cat(paste0("|->",round(100*i/length(idx.x)),"%"), "->")
        if(i>unit & i<length(idx.x)) cat(paste0(round(100*i/length(idx.x)),"%"),"->")
        if(i==length(idx.x))cat("100%|")
    }

}
 cat("\n") 
 rownames(res)=x.names
  if(padjust){
    res[,"dR2.pvalue"]=p.adjust(res[,"dR2.pvalue"], method="BH")
    res[,"dRMSE.pvalue"]=p.adjust(res[,"dRMSE.pvalue"], method="BH")
    res[,"dAccuracy.pvalue"]=p.adjust(res[,"dAccuracy.pvalue"], method="BH")
  }
 
  res=structure(res, class="bprfimp",padjust=padjust)
  return(res)
}

summary.bprfimp=function(bprfimpObj){
  obj=bprfimpObj
  padjust=attr(bprfimpObj, "padjust")
  rr=if(padjust) which(obj[,"dR2.pvalue"]<0.05 & obj[,"dRMSE.pvalue"]<0.05 & obj[,"dAccuracy.pvalue"]<0.05) else {
    which(obj[,"dR2.pvalue"]==min(obj[,"dR2.pvalue"]) & obj[,"dRMSE.pvalue"]==min(obj[,"dRMSE.pvalue"]) & obj[,"dAccuracy.pvalue"]==min(obj[,"dAccuracy.pvalue"]))
  }
  cat("Significant gene blocks:","\n")
  print(obj[rr,,drop=F])
  return(invisible(structure(obj[rr,,drop=F], sig.gene=rr)))
}



