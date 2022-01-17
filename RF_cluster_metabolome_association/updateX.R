accuracy=function(y, yhat){
  return(abs(y-yhat)/abs(y))
}



updateX=function(X, B_pls, cluster, scale.x=NULL){

  
  X0=NULL; if(0 %in% cluster) X0=X[,cluster==0,drop=F]
  x.block=lapply(1:max(cluster), function(cl)X[,cluster==cl,drop=F])
  
    #### scale X according to the X_old that generate B_pls
    if(!is.null(scale.x)){
      x.block.s=lapply(1:length(x.block), function(i){
        x.center=scale.x[[i]]$center;
        x.scale=scale.x[[i]]$scale;
        scale(x.block[[i]], center=x.center, scale=x.scale)})
      x.block=x.block.s  
    }
    
    ##### scale X to one dimension using B_pls
    Xnew=sapply(1:length(x.block), function(i){
      xn=x.block[[i]]
      bb=B_pls[[i]]
      xn %*% bb
    })
    Xnew=matrix(Xnew, nrow(X), max(cluster))
    colnames(Xnew)=paste0("Cluster",1:max(cluster))
    Xnew=cbind(Xnew, X0)
 
  return(Xnew)
}


