
# sPLS 
# Returned values
# UUp
# higher lamda --> higher penalty
sPLS=function(X, Y, lv=1, max.iter=500, tol=1e-06, lamda1=0, lamda2=0, scale=T, verbose=T, quality.check=T){
  ### Data cleaning
  X=as.matrix(X)
  Y=as.matrix(Y)
  lv=min(c(lv, ncol(Y)))
  
  if(is.null(colnames(X))) colnames(X)=paste0("X",1:ncol(X))
  if(is.null(colnames(Y))) colnames(Y)=paste0("Y",1:ncol(Y)) 
  if(is.null(rownames(X))) rownames(X)=paste0("Sample", 1:nrow(X))
  
  if(nrow(X)!=nrow(Y)) stop("X and Y should have equal number of observations")
  H=lv
  keep.x=1:ncol(X)
  keep.y=1:ncol(Y)
  
  if(quality.check){
    if(verbose) cat("Check and Remove X and Y variables with zero or near zero variance ...","\n")
    keep.x=setdiff(1:ncol(X), nearZeroVar(X))
    keep.y=setdiff(1:ncol(Y), nearZeroVar(Y))
    
    X=X[, keep.x, drop=F]
    Y=Y[, keep.y, drop=F]
  }
  
  ### Data transformation
  if(verbose) cat("Scale X and Y variables...","\n")
  Xs=scale(X, center=T, scale=scale)
  Ys=scale(Y, center=T, scale=scale)
  if(scale==F){
    attr(Xs, "scaled:scale")=rep(1, ncol(Xs))
    attr(Ys, "scaled:scale")=rep(1, ncol(Ys))
  } 
  scale.x=list(center=attr(Xs, "scaled:center"),scale=attr(Xs,"scaled:scale"))
  scale.y=list(center=attr(Ys, "scaled:center"),scale=attr(Ys,"scaled:scale"))
  p=ncol(Xs);q=ncol(Ys);n=nrow(Xs)
  
  if(verbose) cat("Initiate sPLS...","\n")
  ## score matrix
  TT=UU=matrix(0, n, H)
  ## weight matrix
  WW=matrix(0, p, H); VV=matrix(0, q, H)
  ## orthogonal loading matrix
  PP=matrix(0, p, H); QQ=matrix(0, q, H)
  ## rotation matrix
  #RR=matrix(0, p, H)
  ## R2Y and R2X
  R2Y=R2X=rep(0, H)

  h=1;Xh=Xs;Yh=Ys
  while(h<=H)
  {
    if(verbose) cat("Search for the weight vectors of X and Y for component ", h, "...","\n")
    wv.h=.sparse.pls(X=Xh, Y=Yh, max.iter=max.iter, tol=tol, lamda1=lamda1, lamda2=lamda2)
    
    if(verbose) cat("   Deflate X and Y","\n")
    ww=wv.h$w.new; vv=wv.h$v.new
    def=.deflate(X=Xh, Y=Yh, ww=ww)
    Xh=def$X.deflate; Yh=def$Y.deflate
    
    if(verbose) cat("   Add score and loading vectors for X and Y to the summary matrix...","\n")
    #X and Y scores for component h
    TT[,h]=def$tt; UU[,h]=def$uu
    #X and Y loadings for component h
    PP[,h]=def$pp; QQ[,h]=def$qq
    #X and Y weight vectors for component h
    WW[,h]=ww; VV[,h]=vv
   
    h=h+1
  }
  
  if(verbose) cat("Model evaluation...", "\n")
  R2=.r2(TT=TT, PP=PP, QQ=QQ, VV=VV, X=Xs, Y=Ys)
  
  ###### calculate B_pls that can be used for prediction #######
  if(verbose) cat("Calculate prediction parameters...", "\n")
  ### Ys =b %*% TT %*% t(QQ) + E, where b is a constant and E is error term
  ### b*t(QQ) = solve(TT, Ys)
  coef = lm(Ys ~ TT)$coefficients
  if(ncol(Ys)==1) BQ=matrix(coef["TT"],1,1) else BQ=matrix(coef[-1,], lv, ncol(Ys))
    
  ## TT = Xs %*% RR
  ## Yhat = b* TT %*% t(QQ) = b* Xs %*% R %*% t(QQ) = Xs %*% R %*% (b*t(QQ))
  ## Let B_pls = R %*% (b*t(QQ))
  ## Yhat = Xs %*% B_pls
  RR = WW %*% solve(t(PP) %*% WW)
  B_pls = RR %*% BQ
  
  colnames(TT)=colnames(UU)=colnames(PP)=colnames(QQ)=colnames(WW)=colnames(VV)=paste0("Comp",1:H)
  rownames(TT)=rownames(UU)=rownames(Xs)
  rownames(PP)=rownames(WW)=rownames(B_pls)=colnames(Xs)
  rownames(QQ)=rownames(VV)=colnames(B_pls)=colnames(Ys)
  
  res.spls=structure(list(TT=TT, UU=UU, PP=PP, QQ=QQ, WW=WW, VV=VV, R2=R2, X=Xs, Y=Ys, keep.x=keep.x, keep.y=keep.y, scale.x=scale.x, scale.y=scale.y, H=H, B_pls=B_pls, lamda1=lamda1, lamda2=lamda2), class="spls")
  if(verbose) cat("Finished", "\n")
  return(invisible(res.spls))
}

vips.spls=function(splsObj){
  WW=splsObj$WW; p=nrow(WW); h=ncol(WW)
  R2Y=splsObj$R2["R2Y",1:h, drop=F]
  total.R2Y=splsObj$R2["R2Y","total"]
  W2=WW^2 ### sum(w2)==1
  vips=sqrt(apply(W2%*%t(R2Y), 1, sum)/total.R2Y)*sqrt(p)
  return(vips)
}

predict.spls=function(splsObj, Xnew, Ynew=NULL){
  Xnew = as.matrix(Xnew)
  if(!is.null(Ynew)) Ynew = as.matrix(Ynew)
  
  x.center = splsObj$scale.x$center
  x.scale = splsObj$scale.x$scale
  
  y.center = splsObj$scale.y$center
  y.scale = splsObj$scale.y$scale
  
  WW = splsObj$WW
  PP = splsObj$PP
  QQ = splsObj$QQ
  B_pls = splsObj$B_pls
  
  Xnew = Xnew[,splsObj$keep.x]
  Xnew = scale(Xnew, center=x.center, scale=x.scale)
  
  Yhat = Xnew %*% B_pls
  if(length(y.scale)==1) y.scale = as.matrix(y.scale) else y.scale = diag(y.scale)
  
  Yhat = Yhat %*% y.scale
  Yhat = sweep(Yhat, 2, STATS = y.center, FUN="+")
  
  if(!is.null(Ynew)){    
    rmse = sqrt(sum((Ynew-Yhat)^2)/nrow(Ynew))
    r2 = mean(sapply(1:ncol(Ynew), function(i)summary(lm(Ynew[,i]~Yhat[,i]))$r.squared)) 
    accuracy =  mean(sapply(1:ncol(Ynew), function(i) abs(Ynew[,i]-Yhat[,i])/abs(Ynew[,i])))
    return(list(Yhat=Yhat, rmse=rmse, r2=r2, accuracy=accuracy))
  }
  
  return(Yhat)
}

.sparse.pls=function(X, Y, p, q, n, max.iter, tol, lamda1, lamda2){
  p=ncol(X); q=ncol(Y); n=nrow(X)
  M = t(X)%*%Y
  svd.M=svd(x=M, nu=1, nv=1)
  w.old=matrix(svd.M$u, p, 1); v.old=matrix(svd.M$v, q, 1)
  
  iter=0
  repeat{
    iter=iter+1
    w.new=matrix(.soft.threshold(M %*% v.old, lamda1), p, 1, dimnames=list(colnames(X), NULL))
    w.new=w.new/sqrt(drop(t(w.new) %*% w.new))
    v.new=matrix(.soft.threshold(t(M) %*% w.new, lamda2), q, 1, dimnames=list(colnames(Y), NULL))
    v.new=v.new/sqrt(drop(t(v.new) %*% v.new))
    
    w.diff=w.old-w.new
    if(drop(t(w.diff)%*%w.diff)<tol | iter>max.iter) return(list(w.new=w.new, v.new=v.new))
    w.old=w.new; v.old=v.new
  } 
}


.soft.threshold=function(input, lamda){
 # input=as.vector(input)
  if(lamda==0) return(input)
  cutoff=quantile(abs(input), lamda)
  output=abs(input)-cutoff
  output[output<=.Machine$double.eps]=0
  return(sign(input)*output)
}  
  
 .deflate=function(X, Y, ww){
   tt = X %*% ww / drop(t(ww) %*% ww) ## X scores
   ## Normalize 
   tt = tt / sqrt(drop(t(tt)%*%tt))
   
   pp = t(X) %*% tt / drop(t(tt)%*%tt) # X loadings
   qq = t(Y) %*% tt / drop(t(tt)%*%tt) # Y loadings 
   uu = Y %*% qq # Y scores calculated based on X weight vector
 
   X.deflate= X - tt %*% t(pp)
   Y.deflate= Y - tt %*% t(qq)
   
   return(list(X.deflate=X.deflate, Y.deflate=Y.deflate, tt=tt, uu=uu, pp=pp, qq=qq))
 } 
    
.r2=function(TT, PP, QQ, VV, X, Y){
   
  H=ncol(TT)
  
  SSX0=sum(X^2); SSY0=sum(Y^2)
  
  R2X=sapply(1:H, function(h){
    X_hat=TT[,h] %*% t(PP[,h]);sum(X_hat^2) / SSX0})
  
  R2Y=sapply(1:H, function(h){
    Y_hat=TT[,h] %*% t(QQ[,h]);sum(Y_hat^2)/SSY0})
  
  R2Y.spar=sapply(1:H, function(h){nonzero.y=which(VV[,h]!=0);SSY1=sum(Y[,nonzero.y]^2);
    Y_spar_hat=TT[,h] %*% t(QQ[nonzero.y,h]); sum(Y_spar_hat^2)/SSY1})
  
  names(R2X)=names(R2Y)=names(R2Y.spar)=paste0("Comp",1:H)
  
  return(cbind(rbind(R2X, R2Y, R2Y.spar), total=c(sum(R2X), sum(R2Y), sum(R2Y.spar)/H)))
}

#sPLS=function(X, Y, lv=1, max.iter=500, tol=1e-06, lamda1=0, lamda2=0, scale=T, verbose=T, quality.check=T){
permute=function(splsObj, B=10000, nThreads=1, p.adjust=T){
  Xs=splsObj$X
  Ys=splsObj$Y
  lamda1=splsObj$lamda1
  lamda2=splsObj$lamda2
  lv=ncol(splsObj$TT)
  scale=T
  if(sum(splsObj$scale.x$scale==1)==ncol(Xs)) scale=F
  vips0=vips.spls(splsObj)
  ss0=splsObj$R2["R2Y","total"]*vips0
  
  permute.one=function(X, Y, lv, lamda1, lamda2, scale){
    X1=X[sample(1:nrow(X), nrow(X)), ,drop=F]
    spb=sPLS(X=X1, Y=Y, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=scale, verbose=F, quality.check=F)
    vipsb=vips.spls(spb)
    ssb=spb$R2["R2Y","total"]*vipsb
    return(ssb>=ss0)
  }
    
  if(nThreads==1){
    pval=NULL
    pb=txtProgressBar(min=0, max=B, char="*", style=3)
    for(b in 1:B){
      pval=cbind(pval, permute.one(X=Xs, Y=Ys, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=scale))
      setTxtProgressBar(pb, b)
    } 
  } else{
     cat("Start parallel permutation...")
     cl=snow::makeCluster(nThreads);doSNOW::registerDoSNOW(cl)
     pb=txtProgressBar(max=B, style=3)
     progress= function(n) setTxtProgressBar(pb, n)
     opts=list(progress = progress)

     pval=foreach(b=1:B, .combine=cbind, .options.snow=opts) %dopar% {
      source("/work/LAS/myn-lab/KChen/SPLS/sPLS.R")
      permute.one(X=Xs, Y=Ys, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=scale)
    }
    snow::stopCluster(cl)
    cat("\n")
  }
  
  count=rowSums(pval) 
  pval=(rowSums(pval)+1)/(1+B)
  if(p.adjust==T){
    padj=p.adjust(pval, method="BH")
    if(sum(padj<0.05)==0) {
	padj=pval
	padj[vips0>1]=p.adjust(padj[vips0>1], method="BH")
     } 
  }
  genes=cbind(importance=vips0, count=count, pval=pval, padj=padj)
  
  return(invisible(genes))
}

tune.spls=function(X, Y, lv=c(1 ,5, 10), scale=T, lamda1=c(0.2, 0.4, 0.6), plot=F, cv.fold=4, cv.repeats=1, nThreads=1, verbose=1,...){
  
  lamda1.spls=rep(lamda1, length(lv))
  ncomp=rep(lv, each=length(lamda1))
  tuned=data.frame(lv=ncomp, lamda1=lamda1.spls)
  
  r2.tune=rmse.tune=accuracy.tune=rep(0, nrow(tuned))
  
  tune.cv=function(X, Y, lv, lamda1, verbose, nThreads, scale){
    X=as.matrix(X)
    Y=as.matrix(Y)
    
    check.x=check.y=1
    
    while(check.x>0 | check.y>0){
      
      folds=cut(sample(1:nrow(X), nrow(X)), cv.fold, labels=F)
      subset.x=lapply(1:cv.fold, function(j)X[folds==j, ,drop=F])
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
      Xtrain=X[-idx, ,drop=F]; Xtest=X[idx, ,drop=F]
      Ytrain=Y[-idx, , drop=F]; Ytest=Y[idx, ,drop=F]
      mod=sPLS(X=Xtrain, Y=Ytrain, lv=lv, lamda1=lamda1, scale=scale, verbose=F, quality.check=F)
      pred=predict(mod, Xnew=Xtest, Ynew=Ytest)
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
      # tune.cv=function(X, Y, lv, lamda1, verbose, nThreads){
      cv.result=tune.cv(X=X, Y=Y, lv=tuned[i,"lv"], lamda1=tuned[i,"lamda1"], verbose=1, nThreads=nThreads, scale=scale)
      r2[r,]=cv.result[,"r2"]
      rmse[r,]=cv.result[,"rmse"]
      accuracy[r,]=cv.result[,"accuracy"]
      cat("\n")
    }
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
  
