  
randomSPLS=function(X, Y, lv=1, lamda1=0, lamda2=0, scale=T,quality.check=T, verbose=T, B=1000, nFeatures=round(sqrt(ncol(X))), nThreads=1){
  X=as.matrix(X)
  Y=as.matrix(Y)
  if(is.null(colnames(X))) colnames(X)=paste0("X",1:ncol(X))
  if(is.null(colnames(Y))) colnames(Y)=paste0("Y",1:ncol(Y)) 
  if(is.null(rownames(X))) rownames(X)=paste0("Sample", 1:nrow(X))
  
  if(nrow(X)!=nrow(Y)) stop("X and Y should have equal number of observations")    
  
  keep.x=1:ncol(X); keep.y=ncol(Y)
  if(quality.check){
    if(verbose) cat("Check and Remove X and Y variables with zero or near zero variance ...","\n")
    keep.x=setdiff(1:ncol(X), nearZeroVar(X))
    keep.y=setdiff(1:ncol(Y), nearZeroVar(Y)) 
    X=X[, keep.x, drop=F]
    Y=Y[, keep.y, drop=F]
  }
  
  randomSPLS.one=function(X, Y, lv, lamda1, lamda2, scale, nFeatures){
    xinB=rep(0, ncol(X))
    xinB[sample(1:ncol(X), nFeatures)]=1
    
    inB=sample(1:nrow(X), nrow(X), replace=T)
    
    while(0 %in% apply(X[inB,xinB==1,drop=F],2,sd) | 0 %in% apply(Y[inB, ,drop=F], 2, sd)){
   	    inB=sample(1:nrow(X), nrow(X), replace=T)
    } 

    ooB=setdiff(1:nrow(X), inB)
    Yhat=matrix(0, nrow(Y), ncol(Y))    
    vips=rep(0, ncol(X))
    Xtrain=X[inB, xinB==1, drop=F]; Xtest=X[ooB, xinB==1, drop=F]
    Ytrain=Y[inB, ,drop=F]; Ytest=Y[ooB, ,drop=F]
    sp=sPLS(X=Xtrain, Y=Ytrain, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=scale, verbose=F, quality.check=F)
    vips[xinB==1]=vips.spls(sp)
    pred=predict(sp, Xtest, Ytest)
    Yhat[ooB,]=pred$Yhat; r2=pred$r2; rmse=pred$rmse; accuracy=pred$accuracy
    
    return(list(tree=sp, Yhat=Yhat, vips=vips, r2=r2, rmse=rmse, accuracy=accuracy, inB=inB, ooB=ooB, xinB=xinB))
  }
  
  if(nThreads==1){
    if(verbose) pb=txtProgressBar(1, B, char="*", style=3)
    spls.tree=list()
    for(b in 1:B){
      spls.tree[[b]]=randomSPLS.one(X=X, Y=Y, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=scale, nFeatures=nFeatures)
      if(verbose){ 
        setTxtProgressBar(pb, b)
        if(b==B)cat("\n")
      }
    }
  } else{
    if(verbose==T) cat("Parallization:\n")
    cl=snow::makeCluster(nThreads);doSNOW::registerDoSNOW(cl) 
     pb=txtProgressBar(max=B, style=3)
     progress= function(n) setTxtProgressBar(pb, n)
     opts=list(progress = progress)
    spls.tree=foreach(b=1:B, .options.snow=opts) %dopar%{
      source("/work/LAS/myn-lab/KChen/RSPLS/sPLS.R")
      randomSPLS.one(X=X, Y=Y, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=scale, nFeatures=nFeatures) 
    }
    snow::stopCluster(cl)
  }
  
  cat("\n")
  spls.tree=do.call(rbind, spls.tree)
  forest=spls.tree[,"tree"]
  vips.tree=do.call(cbind, spls.tree[,"vips"])
  r2.tree=do.call(c, spls.tree[,"r2"])
  rmse.tree=do.call(c, spls.tree[,"rmse"])
  accuracy.tree=do.call(c, spls.tree[,"accuracy"])
  xinB.tree=spls.tree[,"xinB"]
  inB=spls.tree[,"inB"]
  ooB=spls.tree[,"ooB"]
  
  xinB.forest=rowSums(do.call(cbind, xinB.tree))
  xinB.forest[xinB.forest==0]=1
  vips.forest=rowSums(vips.tree)/xinB.forest
  names(vips.forest)=colnames(X)
  xinB=lapply(xinB.tree, function(x)which(x==1))
  
  Yhat=spls.tree[,"Yhat"]
  if(ncol(Y)==1) Yhat.forest=matrix(rowSums(do.call(cbind,Yhat)), nrow(Y), ncol(Y)) else Yhat.forest=Reduce("+", Yhat)
  ooB.fin=rep(0, nrow(Y))
  ooB.table=table(do.call(c, ooB))
  ooB.fin[as.numeric(names(ooB.table))]=ooB.table
  if(ncol(Y)==1) Yhat.forest[ooB.fin>0,]=Yhat.forest[ooB.fin>0,]/ooB.fin[ooB.fin>0] else Yhat.forest[ooB.fin>0,]=sweep(Yhat.forest[ooB.fin>0,], 1, STATS=ooB.fin[ooB.fin>0], FUN="/")
  Y.forest=Y

  #Yhat.forest=apply(Yhat, 1, function(y)mean(y[y!=0]))
  #Y.forest=Y
  #if(sum(is.na(Yhat.forest))>0){
  #  Y.forest=Y[!is.na(Yhat.forest),]
  #  Yhat.forest=Yhat.forest[!is.na(Yhat.forest)]
  #}
  
  r2.forest=mean(sapply(1:ncol(Y.forest[ooB.fin>0,,drop=F]), function(i)summary(lm(Y.forest[ooB.fin>0,i]~Yhat.forest[ooB.fin>0,i]))$r.squared)) 
  rmse.forest = sqrt(sum((Y.forest[ooB.fin>0,] - Yhat.forest[ooB.fin>0,])^2)/nrow(Y.forest[ooB.fin>0,,drop=F]))
  accuracy.forest = mean(sapply(1:ncol(Y.forest[ooB.fin>0,,drop=F]), function(i) abs(Y.forest[ooB.fin>0,i]-Yhat.forest[ooB.fin>0,i])/abs(Y.forest[ooB.fin>0,i])))

  res=structure(list(forest=forest, inB=inB, ooB=ooB, xinB=xinB, X=X, Y=Y, keep.x=keep.x, lamda1=lamda1, lamda2=lamda2, lv=lv, nFeatures=nFeatures,Yhat=Yhat,r2.forest=r2.forest, rmse.forest=rmse.forest, accuracy.forest=accuracy.forest, importance=vips.forest),class="rspls")
     
  return(invisible(res))
}

predict.rspls=function(obj, Xnew, Ynew=NULL){
  keep.x = obj$keep.x
  xinB = obj$xinB
  forest = obj$forest
  
  if(is.null(dim(Ynew))) Ynew=matrix(Ynew, nrow(Xnew), 1)
  
  B = length(forest)
  
  if(ncol(Xnew) > ncol(obj$X)) Xnew = Xnew[,keep.x,drop = F]
  
  Yhat = lapply(1:B, function(b){
    tree = forest[[b]]
    xb = xinB[[b]]
    xn = Xnew[,xb,drop = F]
    predict(splsObj = tree, Xnew = xn)})
   
  #Yhat=rowSums(do.call(cbind, Yhat))/B
  if(ncol(Ynew)==1) Yhat=matrix(rowSums(do.call(cbind, Yhat))/B, nrow(Ynew), ncol(Ynew)) else Yhat=Reduce("+", Yhat)/B

  if(!is.null(Ynew)){
    r2=mean(sapply(1:ncol(Ynew), function(i)summary(lm(Ynew[,i]~Yhat[,i]))$r.squared))
    rmse = sqrt(sum((Ynew - Yhat)^2)/nrow(Ynew))
    accuracy = mean(sapply(1:ncol(Ynew), function(i) abs(Ynew[,i]-Yhat[,i])/abs(Ynew[,i])))
    
    return(list(Yhat=Yhat, r2=r2, rmse=rmse, accuracy=accuracy))
  }
  
  return(list(Yhat=Yhat, na=na))
}

subset.rspls=function(obj, idx.tree){
  subset=list()
  ntree=length(obj$forest)
  for(i in 1:length(obj)){
    if(length(obj[[i]])!=ntree) subset[[i]]=obj[[i]] else subset[[i]]=obj[[i]][idx.tree]
  }
  names(subset)=names(obj)
  return(structure(subset, class="rspls"))
}

importance.one=function(obj, permutation=100, idx.x, nThreads=1){
  ### Find subset trees
  #cat(idx.x)
  xinB=obj$xinB
  idx.tree=which(sapply(xinB, function(x)as.numeric(idx.x %in% x))==1)  
  x.obj=subset(obj, idx.tree=idx.tree)
  X0=x.obj$X
  Y0=x.obj$Y
  Yhat0=obj$Yhat[,idx.tree, drop=F]
  Yhat0=apply(Yhat0, 1, function(y)mean(y[y!=0]))
  
  eval0=c(r2=summary(lm(Y0~Yhat0))$r.squared, rmse=sqrt(sum((Y0-Yhat0)^2)/length(Y0)), accuracy=mean(abs(Y0-Yhat0)/abs(Y0)))
  eval1=matrix(0, permutation, 3, dimnames=list(paste0("P",1:permutation), c("r2","rmse","accuracy")))
  
  perm=function(X0, Y0, x.obj, idx.x){
       idx.pm=sample(1:nrow(X0),nrow(X0),replace=F)
       X1=X0; X1[,idx.x]=X1[idx.pm, idx.x]
       pred=predict.rspls(x.obj, Xnew=X1, Ynew=Y0)
       r2=pred$r2
       rmse=pred$rmse
       accuracy=pred$accuracy
       return(c(r2=r2, rmse=rmse, accuracy=accuracy))  
     }
     
  if(nThreads==1){
	for(p in 1:permutation){
     		eval1[p,]=perm(X0, Y0, x.obj, idx.x)
     		}
	} else{
	cl=snow::makeCluster(nThreads);doSNOW::registerDoSNOW(cl)
    	eval1=foreach(p=1:permutation, .combine=rbind) %dopar% {
  		source("sPLS.R")
      source("randomSPLS.R")
      #source("/work/LAS/myn-lab/KChen/RSPLS/sPLS.R");		
      #source("/work/LAS/myn-lab/KChen/RSPLS/randomSPLS.R");
		perm(X0, Y0, x.obj, idx.x)
	} 
    		snow::stopCluster(cl)
     }
     
  dEval=sweep(eval1, 2, STATS=eval0, FUN="-")
  
  dR2=mean(dEval[,1]); dR2.pvalue=(sum(dEval[,1]>=0)+1)/(1+permutation)
  dRMSE=mean(dEval[,2]); dRMSE.pvalue=(sum(dEval[,2]<=0)+1)/(1+permutation)  
  dAccuracy=mean(dEval[,3]); dAccuracy.pvalue=(sum(dEval[,3]<=0)+1)/(1+permutation) 
     
  return(c(dR2=dR2, dR2.pvalue=dR2.pvalue, dRMSE=dRMSE, dRMSE.pvalue=dRMSE.pvalue,dAccuracy=dAccuracy, dAccuracy.pvalue=dAccuracy.pvalue))
}


importance=function(obj, permutation=100, idx.x, verbose=1, padjust=T, nThreads=1){
  X=obj$X
  x.names=colnames(X)[idx.x]
  res=NULL
  for(i in 1:length(idx.x)){
    res=rbind(res, importance.one(obj, permutation=permutation, idx.x=idx.x[i],nThreads=nThreads))
    unit=ifelse(length(idx.x)>10,round(length(idx.x)/10),1)
    if(verbose>0 & i%%unit==0){
        if(i==unit) cat(paste0("|->",round(100*i/length(idx.x)),"%"), "->")
        if(i>unit & i<length(idx.x)) cat(paste0(round(100*i/length(idx.x)),"%"),"->")
        if(i==length(idx.x))cat("100%|")
    }  
  }
  rownames(res)=x.names
  if(padjust){
    res[,"dR2.pvalue"]=p.adjust(res[,"dR2.pvalue"], method="BH")
    res[,"dRMSE.pvalue"]=p.adjust(res[,"dRMSE.pvalue"], method="BH")
    res[,"dAccuracy.pvalue"]=p.adjust(res[,"dAccuracy.pvalue"], method="BH")
  }
  
  res=structure(res, class="rsplsimp",padjust=padjust)
  return(res)
}

summary.rplsimp=function(obj){
  padjust=attr(obj, "padjust")
  rr=if(padjust) which(obj[,"dR2.pvalue"]<0.05 & obj[,"dRMSE.pvalue"]<0.05 & obj[,"dAccuracy.pvalue"]<0.05) else {
    which(obj[,"dR2.pvalue"]==min(obj[,"dR2.pvalue"]) & obj[,"dRMSE.pvalue"]==min(obj[,"dRMSE.pvalue"]) & obj[,"dAccuracy.pvalue"]==min(obj[,"dAccuracy.pvalue"]))
  }
  cat("Significant gene blocks:","\n")
  print(obj[rr,,drop=F])
  return(invisible(structure(obj[rr,,drop=F], sig.gene=rr)))
}


tune=function(X, Y, lv=1, cluster=rep(1, ncol(X)), lamda1=seq(0, 0.8, 0.2), nFeatures=round(c(ncol(X)/3, ncol(X)/4, sqrt(ncol(X)))),
              plot=T, cv.fold=10, cv.repeats=10, nThreads=1, ...){
  
  l1=rep(lamda1, each=length(nFeatures))
  nF=rep(nFeatures, length(lamda1))
  tuned=data.frame(lamda1=l1, nFeatures=nF)
  r2.tune=rmse.tune=accuracy.tune=rep(0, nrow(tuned))

  for(i in 1:nrow(tuned)){
    if(nThreads==1){
      r2=rmse=accuracy=matrix(0, cv.repeats, cv.fold, dimnames=list(paste0("CV",1:cv.repeats), paste0("Fold",1:cv.fold)))
      for(r in 1:cv.repeats){
        folds=cut(sample(1:nrow(X), nrow(X)), cv.fold, labels=F)
        for(j in 1:cv.fold){
          idx=which(folds==j)
          Xtrain=X[-idx, ,drop=F]; Xtest=X[idx, ,drop=F]
          Ytrain=Y[-idx]; Ytest=Y[idx]
          mod=randomSPLS(X=Xtrain, Y=Ytrain, lv=lv, lamda1=tuned[i,"lamda1"], nFeatures=tuned[i,"nFeatures"], verbose=T, quality.check=F)
          pred=predict(mod, Xtest, Ytest)
          r2[r, j]=pred$r2
          rmse[r, j]=pred$rmse
          accuracy[r, j]=pred$accuracy
        }
      }  
    } else{
      cl=snow::makeCluster(nThreads);doSNOW::registerDoSNOW(cl) 
      pred=foreach(r=1:cv.repeats) %dopar%{
       source("/work/LAS/myn-lab/KChen/RSPLS/sPLS.R")
       source("/work/LAS/myn-lab/KChen/RSPLS/randomSPLS.multi.R")
       folds=cut(sample(1:nrow(X), nrow(X)), cv.fold, labels=F)
       ppred=cbind(r2=rep(0, cv.fold), rmse=rep(0, cv.fold), accuracy=rep(0, cv.fold))
       for(j in 1:cv.fold){
         idx=which(folds==j)
         Xtrain=X[-idx, ,drop=F]; Xtest=X[idx, ,drop=F]
         Ytrain=Y[-idx]; Ytest=Y[idx]
         mod=randomSPLS(X=Xtrain, Y=Ytrain, lv=lv, lamda1=tuned[i,"lamda1"], nFeatures=tuned[i,"nFeatures"], verbose=T, quality.check=F)
         pred=predict(mod, Xtest, Ytest)
         ppred[j,"r2"]=pred$r2
         ppred[j,"rmse"]=pred$rmse
         ppred[j,"accuracy"]=pred$accuracy
       }
     }
     snow::stopCluster(cl)
     r2=t(sapply(pred, function(p)p[,"r2"]))
     rmse=t(sapply(pred, function(p)p[,"rmse"]))
     accuracy=t(sapply(pred, function(p)p[,"accuracy"]))
     colnames(r2)=colnames(rmse)=colnames(accuracy)=paste0("Fold",1:cv.fold)
     rownames(r2)=rownames(rmse)=rownames(accuracy)=paste0("CV",1:cv.repeats)
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
  eval=rbind(r2.tune, rmse.tune, accuracy.tune)
  colnames(eval)=paste(tuned[,"lamda1"], tuned[,"nFeatures"], sep="&")
  
  return(eval) 
}

permute.rspls=function(obj, PP=10000, nThreads=1){
 
  print("Load X and Y") 
  X=obj$X; Y=obj$Y
  lamda1=obj$lamda1; lamda2=obj$lamda2; lv=obj$lv; nFeatures=obj$nFeatures
  B=length(obj$forest)
  
  print("Start Extracting necessary input from the objects")
  imp0=obj$importance
  p=rep(1, length(imp0))
  
  print("Permutation")
  if(nThreads==1){
  #  cat("Progress bar\n") 
    pb=txtProgressBar(min=1, max=PP, style=3, char="*")
    for(pp in 1:PP){
      X1=X[sample(1:nrow(X), nrow(X)),,drop=F]
      Y1=Y
      imp1=randomSPLS(X=X1, Y=Y1, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=T, quality.check=F, verbose=F, nFeatures=nFeatures, B=B)$importance
   #   if((100*pp/PP) %%10 ==0) cat("     Finish estimation: ",round(100*pp/PP,1), "%","\n")
      setTxtProgressBar(pb, pp)
      p=p+as.numeric(imp1>=imp0)
    }
  } else{
    cat("Parallel start", "\n")
     cl=snow::makeCluster(nThreads);doSNOW::registerDoSNOW(cl)
     pb=txtProgressBar(max=PP, style=3)
     progress= function(n) setTxtProgressBar(pb, n)
     opts=list(progress = progress)
     imp1=foreach(pp=1:PP, .combine=cbind, .options.snow=opts) %dopar%{
      source("/work/LAS/myn-lab/KChen/RSPLS/sPLS.R");
      source("/work/LAS/myn-lab/KChen/RSPLS/randomSPLS.multi.R"); 
      
      X1=X[sample(1:nrow(X), nrow(X)),,drop=F];
      Y1=Y;
      randomSPLS(X=X1, Y=Y1, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=T, quality.check=F, verbose=F, nFeatures=nFeatures, B=B)$importance
    }
    
    snow::stopCluster(cl)
    p=cbind(p,apply(imp1, 2, function(ii)as.numeric(ii>=imp0)))
    p=apply(p, 1, sum)
  }
  
  cat("\n")
  names(p)=names(imp0)
  count=p
  p=p/(1+PP)
  padj=p
  padj[imp0>1]=p.adjust(padj[imp0>1], method="BH")
  return(cbind(importance=imp0, count=count, pval=p, padj=padj))
}
