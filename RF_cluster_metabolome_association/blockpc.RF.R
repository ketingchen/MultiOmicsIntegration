#require(randomForest)
#require(ranger)
#source("PLS.R")
#source("updateX.R")

  
block.pc.rf=function(X, Y, lv=1, ntree=1000, mtry=NULL, cluster=rep(1, ncol(X)), lamda=0, quality.check=T,verbose=1, package=c("ranger","randomForest")){
  if(length(package)==2) package="ranger"
  X=as.matrix(X)
  if(is.null(colnames(X))) colnames(X)=paste0("X", 1:ncol(X))
  
  keep.x=1:ncol(X)
  if(quality.check){
    if(verbose>0) cat("Remove X varialbes with zero or near zero variance...","\n")
    keep.x=setdiff(1:ncol(X), nearZeroVar(X))
    X=X[,keep.x,drop=F]
    cluster=cluster[keep.x]
  }
  
  if(verbose>0) cat("Generate X block variables...","\n")
  cluster.size=max(cluster)
  X0=NULL;if(0 %in% cluster) X0=X[,cluster==0,drop=F]
  x.block=lapply(1:cluster.size, function(cl)X[,cluster==cl,drop=F])
  sp=lapply(1:cluster.size, function(i)sPCA(X=x.block[[i]], lv=lv, scale=T, lamda=lamda, quality.check=F))
  sp=do.call(rbind, sp)
  B_pca=sp[,"PP"]
  scale.x=sp[,"scale.x"]
  R2X = do.call(cbind, sp[,"R2X"])
  Xnew=do.call(cbind,sp[,"TT"])
  if(lv>1) colnames(Xnew)=paste(rep(paste0("Cluster",1:cluster.size), each=lv),1:lv, sep="_") else colnames(Xnew)=paste0("Cluster",1:cluster.size)
  colnames(R2X)=colnames(Xnew)
  Xnew=cbind(Xnew, X0)
  
  if(verbose>0) cat("Build blockPLS Trees...","\n")
  
  bpTree=list()
  for(i in 1:ntree){
    unit=round(ntree/10)
    if(verbose>0 & i%%unit==0){
      if(i%%unit==0){
        if(i==unit) cat("|->10% ->")
        if(i>unit & i<ntree) cat(paste0(100*i/ntree,"%"),"->")
	if(i==ntree)cat("100%|")	
      } 
    }	
    bpTree[[i]]=block.pc.tree(X=Xnew,Y=Y, mtry=mtry, package=package)
  }
  
  if(verbose>0) {
    cat("\n")
    cat("Prepare data summary...","\n")
  }
  
  bpTree=do.call(rbind, bpTree)
  forest=bpTree[,"tree"]
  inB=bpTree[,"inB"]
  ooB=bpTree[,"ooB"]
  xinB=bpTree[,"xinB"]
  Xtest=bpTree[,"Xtest"]
  Yhat=do.call(cbind, bpTree[,"Yhat"])
  r2.tree=do.call(c,bpTree[,"r2"])
  rmse.tree=do.call(c,bpTree[,"rmse"])
  
  Yhat.fin=apply(Yhat, 1, function(yy)mean(yy[yy!=0]))
  Y.fin=Y
  if(sum(is.na(Yhat.fin))>0) {
    Y.fin=Y.fin[!is.na(Yhat.fin)]
    Yhat.fin=Yhat.fin[!is.na(Yhat.fin)]
  }

  rmse.forest=sqrt(sum((Y.fin-Yhat.fin)^2)/length(Y.fin))
  accuracy=abs(Y.fin-Yhat.fin)/abs(Y.fin)
  r2.forest=summary(lm(Y.fin~Yhat.fin))$r.squared
  names(forest)=names(inB)=names(ooB)=names(xinB)=names(Xtest)=colnames(Yhat)=names(r2.tree)=names(rmse.tree)=paste0("Tree",1:ntree)
  
  if(package=="ranger"){
    imp=ranger(Y~., data=data.frame(Xnew, Y=Y), num.tree=ntree, importance="permutation")$variable.importance 
    imp=matrix(imp,dimnames=list(colnames(Xnew)))  
  } else{
    imp=randomForest(Y~., data=data.frame(Xnew, Y=Y), ntree=ntree, importance=T)$importance
  } 
  
  if(verbose>0) cat("Job finished.","\n")
  res=structure(list(forest=forest, inB=inB,ooB=ooB, xinB=xinB, B_pca=B_pca, cluster=cluster, X.block=Xnew, Y=Y, keep.x=keep.x, 
                     Xtest=Xtest, Yhat=Yhat, Yhat.fin=Yhat.fin, pca.R2X=R2X, r2.tree=r2.tree, rmse.tree=rmse.tree,
                     rmse.forest=rmse.forest, r2.forest=r2.forest, accuracy=accuracy, scale.x=scale.x, importance=imp,
                     package=package),class="bprf")
  rm(bpTree)
  return(invisible(res))
}

block.pc.tree=function(X, Y, mtry,package){
  if(is.null(mtry)) mtry=round(sqrt(ncol(X)))
  inB=sample(1:nrow(X), nrow(X), replace=T)
  ooB=setdiff(1:nrow(X), inB)
  xinB=sample(1:ncol(X),mtry,replace=F)
  
  Xtrain=X[inB, xinB, drop=F]
  Xtest=X[ooB, xinB, drop=F]
  Ytrain=Y[inB]
  Ytest=Y[ooB]
  Yhat=rep(0, nrow(X))
  
  rownames(Xtrain)=1:nrow(Xtrain)
  train.data=data.frame(Xtrain, Y=Ytrain)
  if(package=="ranger"){
    tree=ranger(Y~., data=train.data, mtry=mtry, num.trees=1, write.forest=T)
    pred=predict(tree, data=data.frame(Xtest))$predictions    
  } else {
    tree=randomForest(Y~., data=train.data, mtry=mtry, ntree=1)
    pred=predict(tree, newdata=Xtest)
  }
  
  Yhat[ooB]=pred
  r2=summary(lm(Ytest~pred))$r.squared
  rmse=sqrt(sum((Ytest-pred)^2)/length(Ytest))
  return(list(tree=tree, inB=inB, ooB=ooB, xinB=xinB, Xtest=Xtest, Yhat=Yhat, r2=r2, rmse=rmse))
}



predict.bprf=function(bprfObj, Xnew, Ynew=NULL){
 obj=bprfObj
 forest=obj$forest
 B_pca=obj$B_pca
 keep.x=obj$keep.x
 scale.x=obj$scale.x
 cluster=obj$cluster
 package=obj$package
 
 Xnew=Xnew[,keep.x,drop=F]
 X.block.new=updateX(Xnew, cluster=cluster, B_pls=B_pca, scale.x=scale.x) 
 
 Yhat=if(package=="ranger") {
   sapply(forest, function(tree) predict(tree, data=data.frame(X.block.new))$predictions)
  } else sapply(forest, function(tree) predict(tree, newdata=X.block.new))
 
 Yhat.fin=apply(Yhat,1,function(y)mean(y[y!=0]))
 if(!is.null(Ynew)){
   Ynew.fin=Ynew
   if(sum(is.na(Yhat.fin))>0){
     Ynew.fin=Ynew[!is.na(Yhat.fin)]
     Yhat.fin=Yhat.fin[!is.na(Yhat.fin)]
   }
   rmse=sqrt(sum((Ynew.fin-Yhat.fin)^2)/length(Ynew))
   r2=summary(lm(Ynew.fin~Yhat.fin))$r.squared
   accuracy=mean(abs(Ynew.fin-Yhat.fin)/abs(Ynew.fin))
   return(list(Yhat=Yhat.fin, rmse=rmse,r2=r2, accuracy=accuracy))
 }
 return(Yhat.fin)
}


tuneLamda=function(X, Y, lv=1, cluster=rep(1, ncol(X)), lamda=seq(0, 0.8, 0.2), plot=T,cv.fold=10, cv.repeats=10, verbose=1, package=c("ranger","randomForest"), quality.check=F){
  
  if(length(package)==2) package="randomForest"
  
  r2.lam=rmse.lam=accuracy.lam=rep(0, length(lamda))
  
  
  if(verbose>0) cat("Generate X block variables...","\n")
  if(quality.check){
    if(verbose>0) cat("Remove X varialbes with zero or near zero variance...","\n")
    keep.x=setdiff(1:ncol(X), nearZeroVar(X))
    X=X[,keep.x,drop=F]
    cluster=cluster[keep.x]
  }
  
  for(l in 1:length(lamda)){
    lam=lamda[l]
    if(verbose>0) cat("Tune Lamda:", lam, "\n")
    r2=rmse=accuracy=matrix(0, cv.repeats, cv.fold, dimnames=list(paste0("CV",1:cv.repeats), paste0("Fold",1:cv.fold)))
    for(i in 1:cv.repeats){
        folds=cut(sample(1:nrow(X), nrow(X)), cv.fold, labels=F)
        if(verbose>0) cat("Cross-validation: repeat", i, "\n")
        pb=txtProgressBar(min = 1, max = cv.fold, style=3, char="*")
        for(j in 1:cv.fold){
          idx=which(folds==j)
          Xtrain=X[-idx, ,drop=F]; Xtest=X[idx, ,drop=F]
          Ytrain=Y[-idx]; Ytest=Y[idx]
          mod=block.pc.rf(X=Xtrain, Y=Ytrain, lv=lv, cluster=cluster, lamda=lam, verbose=0, package=package, quality.check=F)
          pred=predict.bprf(mod, Xtest, Ytest)
          r2[i, j]=pred$r2
          rmse[i, j]=pred$rmse
          accuracy[i, j]=pred$accuracy
          setTxtProgressBar(pb, j)
        }
        if(verbose>0) cat("\n")
    }
    r2.lam[l]=mean(r2) 
    rmse.lam[l]=mean(rmse) 
    accuracy.lam[l]=mean(accuracy)
  }
  
  if(plot){
    par(mfrow=c(1, 3))
    plot(r2.lam, type="o", pch=16); mtext("R2")
    plot(rmse.lam, type="o", pch=16); mtext("RMSE")
    plot(accuracy.lam, type="o", pch=16); mtext("Accuracy")
    par(mfrow=c(1, 1))
  }
  eval=rbind(r2.lam, rmse.lam, accuracy.lam)
  colnames(eval)=lamda
  rownames(eval)=c("r2","rmse","accuracy")
  
  return(eval) 
}
