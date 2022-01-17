
DataToBlocks=function(X, cluster=rep(1, ncol(X)), lamda=rep(0, max(cluster)), lv=1){
  cluster.size=max(cluster)
  if(cluster.size==0) return(list(Xnew=X))
  
  X0=NULL;if(0 %in% cluster) X0=X[,cluster==0,drop=F]
  x.block=lapply(1:cluster.size, function(cl)X[,cluster==cl,drop=F])
  sp=lapply(1:cluster.size, function(i)sPCA(X=x.block[[i]], lv=lv, scale=T, lamda=lamda[i], quality.check=F))
  sp=do.call(rbind, sp)
  B_pca=sp[,"PP"]
  scale.x=sp[,"scale.x"]
  R2X = do.call(cbind, sp[,"R2X"])
  Xnew =do.call(cbind,sp[,"TT"])
  if(lv>1) colnames(Xnew)=paste(rep(paste0("Cluster",1:cluster.size), each=lv),1:lv, sep="_") else colnames(Xnew)=paste0("Cluster",1:cluster.size)
  colnames(R2X)=colnames(Xnew)
  Xnew=cbind(Xnew, X0)
  return(list(R2X=R2X, Xnew=Xnew, B_pca=B_pca))
}

tuneBlk=function(X, Y, lv=1, cluster=rep(1, ncol(X)), lamda=seq(0, 0.8, 0.2), lamda.rspls=0, nThreads=1, method=c("blockpc.rf","rspls"), verbose=1,...){
  
 # nF=rep(nFeatures, length(lamda1))
 # tuned=data.frame(lamda1=l1, nFeatures=nF)
   r2.tune=rmse.tune=accuracy.tune=matrix(0, nrow=max(cluster), ncol=length(lamda))
    for(cl in 1:max(cluster)){
      unit=round(max(cluster)/10)
      if(verbose>0 & cl%%unit==0){
        if(cl==unit) cat("|->10% ->")
        if(cl>unit & cl<max(cluster)) cat(paste0(round(100*cl/max(cluster)),"%"),"->")
        if(cl==max(cluster)) cat("100%|")
      } 
       for(l in 1:length(lamda)){
         lamda.pca=rep(0, max(cluster))
         lamda.pca[cl]=lamda[l]
         if(method=="rspls"){
           Xnew=DataToBlocks(X, cluster=cluster, lamda=lamda.pca, lv=lv)$Xnew
           mod=randomSPLS(X=Xnew, Y=Y, lv=lv, lamda1=lamda.rspls, quality.check=F, verbose=F)
           r2.tune[cl, l]=mod$r2.forest
           rmse.tune[cl, l]=mod$rmse.forest
           accuracy.tune[cl, l]=mod$accuracy.forest
         } else{
           mod=block.pc.rf(X=X, Y=Y, lv=1, cluster=cluster, quality.check=F, lamda=lamda.pca, package="randomForest", verbose=0)
           r2.tune[cl, l]=mod$r2.forest
           rmse.tune[cl, l]=mod$rmse.forest
           accuracy.tune[cl, l]=mean(mod$accuracy)
         }
       }
    }
    colnames(r2.tune)=colnames(rmse.tune)=colnames(accuracy.tune)=lamda
    rownames(r2.tune)=rownames(rmse.tune)=rownames(accuracy.tune)=paste("Cluster",1:max(cluster), sep="")
    
    return(list(r2=r2.tune, rmse=rmse.tune, accuracy=accuracy.tune))
}

tuneBlk.cv.rspls=function(X, Y, lv=c(1 ,5, 10),cluster=rep(0, ncol(X)), lamda.pca=seq(0.2, 0.4, 0.6), lamda.rspls=c(0.2, 0.4, 0.6), nFeatures=1,plot=F, cv.fold=4, cv.repeats=1, nThreads=1, verbose=1,...){
  
  Xnew=DataToBlocks(X, cluster=cluster)$Xnew
  nFeatures123=c(nFeatures1=round(sqrt(ncol(Xnew))), nFeatures2=round(ncol(Xnew)/2),nFeatures3=round(ncol(Xnew)/3))
  nFeatures=nFeatures123[nFeatures] 
  
  nF=rep(nFeatures, length(lamda.rspls)*length(lamda.pca)*length(nFeatures))
  l.pca=rep(rep(lamda.pca, each=length(nFeatures)), length(lamda.rspls))
  l.rspls=rep(lamda.rspls, each=length(lamda.pca)*length(nFeatures))
  
  ncomp=rep(lv, each=length(l.rspls))
  nF=rep(nF, length(lv))
  l.pca=rep(l.pca, length(lv))
  l.rspls=rep(l.rspls, length(lv))

  tuned=data.frame(lv=ncomp, lamda.rspls=l.rspls, lamda.pca=l.pca, nFeatures=nF)
 
  r2.tune=rmse.tune=accuracy.tune=rep(0, nrow(tuned))
  
  tune.cv=function(X, Y, lv, cluster, lamda.pca, lamda.rspls, nFeatures, verbose, nThreads){
    X=as.matrix(X)
    Y=as.matrix(Y)
    Xnew=DataToBlocks(X, cluster=cluster, lamda=rep(lamda.pca, max(cluster)), lv=1)$Xnew
   
    check.x=check.y=1
    
    while(check.x>0 | check.y>0){

	folds=cut(sample(1:nrow(Xnew), nrow(Xnew)), cv.fold, labels=F)
   	subset.x=lapply(1:cv.fold, function(j)Xnew[folds==j, ,drop=F])
    	subset.y=lapply(1:cv.fold, function(j)Y[folds==j, ,drop=F])
          
      	check.x=sum(sapply(subset.x, function(xx) sum(apply(xx, 2, sd)==0)>0))
        check.y=sum(sapply(subset.y, function(yy) sum(apply(yy, 2, sd)==0)>0))
    }

    eval=NULL
    for(j in 1:cv.fold){
      if(verbose!=0) {
        if(j==1) cat("  Fold")
        cat("_",j)
      }
      idx=which(folds==j)
      Xtrain=Xnew[-idx, ,drop=F]; Xtest=Xnew[idx, ,drop=F]
      Ytrain=Y[-idx, , drop=F]; Ytest=Y[idx, ,drop=F]
      mod=randomSPLS(X=Xtrain, Y=Ytrain, lv=lv, lamda1=lamda.rspls, nFeatures=nFeatures, verbose=F, quality.check=F, nThreads=nThreads)
      pred=predict(mod, Xtest, Ytest)
      eval=rbind(eval, c(r2=pred$r2, rmse=pred$rmse, accuracy=pred$accuracy))
    }
    return(eval)
  }
  
  cat("Total tuning sets:",nrow(tuned),"\n")
  for(i in 1:nrow(tuned)){
      cat("Tuning set:", i,"\n")
      r2=rmse=accuracy=matrix(0, cv.repeats, cv.fold, dimnames=list(paste0("CV",1:cv.repeats), paste0("Fold",1:cv.fold)))
      for(r in 1:cv.repeats){
        cat("   CV repeats:", r)
        cv.result=tune.cv(X=X, Y=Y, lv=tuned[i,"lv"], cluster=cluster, 
                          lamda.pca=tuned[i,"lamda.pca"], lamda.rspls=tuned[i,"lamda.rspls"], 
                          nFeatures=tuned[i,"nFeatures"], verbose=1, nThreads=nThreads)
	r2[r,]=cv.result[,"r2"]
        #print(r2)
	rmse[r,]=cv.result[,"rmse"]
        accuracy[r,]=cv.result[,"accuracy"]
        cat("\n")
      }
   # print(c(r2=mean(r2),rmse=mean(rmse),accuracy=mean(accuracy)))
   # cat("\n")
    r2.tune[i]=mean(r2) 
    rmse.tune[i]=mean(rmse) 
    accuracy.tune[i]=mean(accuracy)
  }
  
  if(plot){
    par(mfrow=c(1, 3))
    plot(r2.tune, type="o", pch=16); mtext("R2")
    plot(rmse.tune, type="o", pch=16); mtext("RMSE")
    plot(accuracy.tune, type="o", pch=16); mtext("Accuracy")
    par(mfrow=c(1, 1))
  }
  cat("\n")
  
  eval=cbind(tuned, r2.tune, rmse.tune, accuracy.tune)
  return(eval)   
}


tuneBlk.bs.rspls=function(X, Y, lv=c(1, 5, 10), cluster=rep(0, ncol(X)), lamda.pca=c(0, 0.4, 0.6), lamda.rspls=c(0, 0.4, 0.6), nFeatures=ncol(X),plot=F, cv.fold=4, cv.repeats=1, nThreads=1, verbose=1,...){
  
  if(length(nFeatures==1)){
    Xnew=DataToBlocks(X, cluster=cluster)$Xnew
    nFeatures=round(sqrt(ncol(Xnew)))
  }
  
  nF=rep(nFeatures, length(lamda.rspls)*length(lamda.pca)*length(nFeatures))
  l.pca=rep(rep(lamda.pca, each=length(nFeatures)), length(lamda.rspls))
  l.rspls=rep(lamda.rspls, each=length(lamda.pca)*length(nFeatures))
  
  ncomp=rep(lv, each=length(l.rspls))
  nF=rep(nF, length(lv))
  l.pca=rep(l.pca, length(lv))
  l.rspls=rep(l.rspls, length(lv))
  
  tuned=data.frame(lv=ncomp, lamda.rspls=l.rspls, lamda.pca=l.pca, nFeatures=nF)
  
  r2.tune=rmse.tune=accuracy.tune=rep(0, nrow(tuned))
  
  for(i in 1:nrow(tuned)){
    if(i==1) cat("Total tuning sets:", nrow(tuned),)
    cat("_",i)
    if(i==nrow(tuned))
      Xnew=DataToBlocks(X, cluster=cluster, lamda=rep(tune[i,"lamda.pca"], max(cluster)), lv=1)$Xnew
    mod=randomSPLS(X=Xnew, Y=Y, lv=tuned[i,"lv"], lamda1=tuned[i,"lamda.rspls"], nFeatures=tuned[i,"nFeatures"], verbose=F, quality.check=F)
    r2.tune[i]=mod$r2.forest
    rmse.tune[i]=mod$rmse.forest
    accuracy.tune[i]=mod$accuracy.forest
  }
  
  if(plot){
    par(mfrow=c(1, 3))
    plot(r2.tune, type="o", pch=16); mtext("R2")
    plot(rmse.tune, type="o", pch=16); mtext("RMSE")
    plot(accuracy.tune, type="o", pch=16); mtext("Accuracy")
    par(mfrow=c(1, 1))
  }
  cat("\n")
  
  eval=cbind(tuned, r2.tune, rmse.tune, accuracy.tune)
  return(eval)   
}
