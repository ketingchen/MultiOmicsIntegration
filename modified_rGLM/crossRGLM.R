require(randomGLM)

#require(doParallel)

.run.rglm=function(X,Y,testIndx=NULL, nThreads,verbose=0)
{
  xtest=X[testIndx,]
  if(length(testIndx)==1) xtest=rbind(X[testIndx,], rep(0, dim(X)[2]))
  ytest=Y[testIndx]
  if(is.null(testIndx)) xtest=ytest=NULL
  if(!is.null(testIndx))
  { 
    x=X[-testIndx,]
    y=Y[-testIndx]
  }
  else
  { 
    x=X; y=Y
  }

  N=dim(X)[2]
  if(verbose>0) cat("Start RGLM...")
  RGLM=randomGLM(x, y, xtest, classify=F, nFeaturesInBag=N*0.2, keepModels=T, nBags=100, nThreads=nThreads)
  if(verbose>0) cat("Finish\n")
  varImp=as.vector(RGLM$timesSelectedByForwardRegression)
  names(varImp)=colnames(x)
  y_hat=RGLM$predictedOOB.response
  estimate=data.frame(y, y_hat)
  if(!is.null(testIndx))
  	{
      predictedTest= RGLM$predictedTest.response
  	  if(length(testIndx)==1) predictedTest=predictedTest[1]
	   # Ifelse statement for k-fold cross validation
      press=ifelse(length(testIndx)>1, sum((ytest-predictedTest)^2)/sum(ytest^2), (ytest-predictedTest)^2/(ytest^2))
      return(list(varImp=varImp, press=press, estimate=estimate))
    }
  rm(RGLM)
  return(list(varImp=varImp, estimate=estimate))
}

# compute the randomGLM results with entire set of data (without subsetting for training and validating sets)
# B: number of permutation runs
whole.rglm=function(X, Y, B=0, numCores=1)
{
  if(B<0) stop("B has to be >= 0.")
  if(B==0) return(.run.rglm(X=X, Y=Y, testIndx=NULL, nThreads=numCores))
  else{
    X0=X;Y0=Y1=Y
    M=dim(X0)[1]
    X1=lapply(1:B, function(b)X0[sample(1:M, M),])
    cat("RGLM with original data\n")
    r0=.run.rglm(X=X0, Y=Y0, testIndx=NULL, nThreads=numCores, verbose=1)
    genes0=r0$varImp;estimate0=r0$estimate
    rm(r0)
  #  counter <- 0
    cat("RGLM with permuted data\n")
    genes1=list()
    for(b  in 1:B){
	cat(b,"_")
   # counter = counter + 1;
    #if(counter %in% seq(0, B, B/10)) print(paste(100*counter/B, "% has been processed"))
    genes1[[b]]=.run.rglm(X=x, Y=Y1, testIndx=NULL, nThreads=numCores, verbose=0)$varImp	   }
    genes1=do.call(cbind, genes1)
    return(list(genes0=genes0, genes1=genes1, estimate0=estimate0))
  }
}


cross.rglm=function(X, Y, subset=c("LeaveOneOut","K-fold"),k=0, numCores=1)
{
 
  #  press=NULL
  #	corOOB=corTest=NULL
  #	varImp=NULL
  #	N=length(vip)  
  if(length(subset)==2)subset=subset[1]
  if(subset=="K-fold"&k==0) stop("k has to be a positive integer for K-fold validation")
  if(subset=="LeaveOneOut"&k>0) stop("No need specify k for Leave-One-Out validation")
  if(!subset %in% c("LeaveOneOut","K-fold")) stop("Subset method has to be Leave-One-Out or K-fold")
  if(dim(X)[1]!=length(Y)) stop("X and Y should have equal number of observations")
  
  if(sum(is.na(Y))>0)
  {
    cat("Warning: samples without response data are removed.", "\n")
    rmv=which(is.na(Y))
    X=X[-rmv,];Y=Y[-rmv]
  }
  
  size=dim(X)[1];seq.new=seq=1:size; brks=size
  if(subset=="K-fold")
  {
    seq.new=sample(1:size, size, replace=F); brks=k
  }
  folds=cut(seq.new, breaks=brks, labels=F)
  
  # function for a single run for cross validation rglm
 # run.rglm=function(X,Y,testIndx)
#  {
#    xtest=X[testIndx,]
#    if(length(testIndx)==1) xtest=rbind(X[testIndx,], rep(0, dim(X)[2]))
#    ytest=Y[testIndx]
#    x=X[-testIndx,]
#    y=Y[-testIndx]
#    N=dim(X)[2]
#    RGLM=randomGLM(x, y, xtest, classify=F, nFeaturesInBag=N*0.2, keepModels=T, nBags=100, nThreads=1)
#    predictedTest= RGLM$predictedTest.response
#    if(length(testIndx)==1) predictedTest=predictedTest[1]
    # Ifelse statement for k-fold cross validation
#    press=ifelse(length(testIndx)>1, sum((ytest-predictedTest)^2)/sum(ytest^2), (ytest-predictedTest)^2/(ytest^2))
#    varImp=as.vector(RGLM$timesSelectedByForwardRegression)
#    names(varImp)=colnames(x)
#    return(list(varImp=varImp, press=press))
#  }
 
  comp=list()
  
  for(i in 1:brks)
  {
    comp[[i]]=.run.rglm(X=X, Y=Y, testIndx=seq[folds==i], nThreads=numCores)
  }

    # using parallel argument to determine whether start parallelization, usually parallel is T unless when the corss.rglm function is called by mc.simulation
#  if(parallel) 
#    {
#      cl=parallel::makeCluster(numCores); doParallel::registerDoParallel(cl)
#      comp=foreach(i=1:brks, .packages="randomGLM") %dopar% run.rglm(X=X,Y=Y,testIndx=seq[folds==i])
#      stopImplicitCluster(); parallel::stopCluster(cl)
#    }
  
#  else
#      comp=foreach(i=1:brks, .packages="randomGLM") %do% run.rglm(X=X,Y=Y,testIndx=seq[folds==i])
  
  comp=do.call(rbind, comp)
  press=do.call(c, comp[,"press"]) # do.call works nicely with the lists
  gene=do.call(cbind, comp[,"varImp"])
 
  if(subset=="K-fold")
    colnames(gene)=paste("Fold",1:k,sep=" ")
  else
    colnames(gene)=paste("LOO",1:dim(X)[1],sep=" ")
  return(list(press=press, gene=gene))
}

 
 
 permute=function(X, Y, B=100, indx=0, subset=c("LeaveOneOut", "K-fold"), k=0, numCores=1){
   
   Y1=Y
   if(indx==0) 
   {
     X1=X[sample(1:dim(X)[1], dim(X)[1]),]
     comp1=cross.rglm(X=X1, Y=Y1, subset=subset, k=k, numCores=numCores)
     return(comp1)
   }
   #if(indx>0)
   #{
  #   X1=lapply(1:B, function(b)X[sample(1:dim(X)[1], dim(X)[1]),])
     
   #}   
  #   X1=X[sample(1:dim(X)[1], dim(X)[1]),]
  # else
  #  X1=cbind(X[, -indx], X[sample(1:dim(X)[1], dim(X)[1]), indx])   
 }
 
 
#monte carlo simulation
mc.simulation=function(X, Y, B=100, subset=c("LeaveOneOut", "K-fold"), k=0, numCores=1, parallel=T)
{
  if(length(subset)==2)subset=subset[1]
  if(subset=="K-fold"&k==0) stop("k has to be a positive integer for K-fold validation")
  if(subset=="LeaveOneOut"&k>0) stop("No need specify k for Leave-One-Out validation")
  if(!subset %in% c("LeaveOneOut","K-fold")) stop("Subset method has to be Leave-One-Out or K-fold")
  
  if(sum(is.na(Y))>0)
  {
    cat("Warning: samples without response data are removed.", "\n")
    rmv=which(is.na(Y))
    X=X[-rmv,];Y=Y[-rmv]
  }
  
  size=dim(X)[1];seq.new=seq=1:size; brks=size
  if(subset=="K-fold")
  {
    seq.new=sample(1:size, size, replace=F); brks=k
  }
  folds=cut(seq.new, breaks=brks, labels=F)
  
  xt=X;yt=Y
  rng.x=apply(xt, 2, range);rng.y=range(yt)
  xr=mclapply(1:B, function(b)apply(rng.x, 2, function(r)runif(size,min=r[1],max=r[2])), mc.cores=numCores)
  yr=mclapply(1:B, function(b)runif(size,min=rng.y[1],max=rng.y[2]), mc.cores=numCores)
  #yr=mclapply(1:B, function(b)yt)
  
  run.rglm=function(X,Y,testIndx)
  {
    xtest=X[testIndx,]
    if(length(testIndx)==1) xtest=rbind(X[testIndx,], rep(0, dim(X)[2]))
    ytest=Y[testIndx]
    x=X[-testIndx,]
    y=Y[-testIndx]
    N=dim(X)[2]
    RGLM=randomGLM(x, y, xtest, classify=F, nFeaturesInBag=N*0.2, keepModels=T, nBags=100, nThreads=1)
    predictedTest= RGLM$predictedTest.response
    if(length(testIndx)==1) predictedTest=predictedTest[1]
    # Ifelse statement for k-fold cross validation
    press=ifelse(length(testIndx)>1, sum((ytest-predictedTest)^2)/sum(ytest^2), (ytest-predictedTest)^2/(ytest^2))
    varImp=as.vector(RGLM$timesSelectedByForwardRegression)
    names(varImp)=colnames(x)
    return(list(varImp=varImp, press=press))
  }
  
  if(parallel) 
  {
    cl=parallel::makeCluster(numCores); doParallel::registerDoParallel(cl)
    ref=foreach(b=1:B, .packages="randomGLM") %:% foreach(i=1:brks, .packages="randomGLM") %dopar% run.rglm(X=X,Y=Y,testIndx=seq[folds==i])
    stopImplicitCluster(); parallel::stopCluster(cl)
  }
  
  else
    ref=foreach(b=1:B, .packages="randomGLM") %:% foreach(i=1:brks, .packages="randomGLM") %do% run.rglm(X=X,Y=Y,testIndx=seq[folds==i])
  
  press.r=mclapply(ref, function(r)do.call(c,lapply(r, function(rr)rr$press)), mc.cores=numCores)
  press.r=sapply(press.r, mean)
  gene.r=mclapply(ref, function(r)lapply(r, function(rr)rr$varImp), mc.cores=numCores)
  gene.r=do.call(rbind, gene.r)
  gene.tmp=apply(gene.r, 1, function(g)rowSums(t(do.call(rbind, g))))
  gene.r=gene.tmp; rm(gene.tmp)
  colnames(gene.r)=paste0("Ref",1:B);rownames(gene.r)=colnames(X)
  return(list(press.ref=press.r, gene.ref=gene.r))
}
