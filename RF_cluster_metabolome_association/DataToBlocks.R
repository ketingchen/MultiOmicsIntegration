source("sPCA.R")
DataToBlocks=function(X, cluster=rep(1, ncol(X)), lamda=rep(0, max(cluster)), lv=1){
  cluster.size=max(cluster)
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

tuneBlk=function(X, Y, lv=1, cluster=rep(1, ncol(X)), lamda=seq(0, 0.8, 0.2), lamda.rspls=0,
                 plot=T, cv.fold=10, cv.repeats=10, nThreads=1, method=c("blockpc.rf","rspls"), verbose=1,...){
  
 # l1=rep(lamda1, each=length(nFeatures))
 # nF=rep(nFeatures, length(lamda1))
 # tuned=data.frame(lamda1=l1, nFeatures=nF)
   r2.tune=rmse.tune=accuracy.tune=matrix(0, nrow=max(cluster), ncol=length(lamda))
    for(cl in 1:max(cluster)){
      unit=round(max(cluster)/10)
      if(verbose>0 & cl%%unit==0){
        if(cl==unit) cat("|->10% ->")
        if(cl>unit & cl<max(cluster)) cat(paste0(100*cl/max(cluster),"%"),"->")
        if(cl==max(cluster))cat("100%|")
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
           mod=block.pc.rf(X=X, Y=Y, lv=1, cluster=cluster, quality.check=F, lamda=lamda.pca, package="randomForest")
           r2.tune[cl, l]=mod$r2.forest
           rmse.tune[cl, l]=mod$rmse.forest
           accuracy.tune[cl, l]=mean(mod$accuracy)
         }
       }
    }
}


