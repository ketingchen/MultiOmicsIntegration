sPCA=function(X, lv=1, max.iter=500, tol=1e-06, lamda=0, scale=T, verbose=T, quality.check=T){
  X = as.matrix(X)
  if(is.null(colnames(X))) colnames(X) = paste0("X",1:ncol(X))
  if(is.null(rownames(X))) rownames(X) = paste0("Sample", 1:nrow(X))
  H = lv
  keep.x=1:ncol(X)
  if(quality.check){
    keep.x = setdiff(1:ncol(X), nearZeroVar(X))
    X=X[, keep.x, drop=F]
  }
  Xs = scale(X, center=T, scale=scale)
  scale.x = list(center=attr(Xs, "scaled:center"),scale=attr(Xs,"scaled:scale"))
  p = ncol(Xs); n = nrow(Xs)
  
  ## score matrix
  TT = matrix(0, n, H, dimnames=list(rownames(Xs), paste0("Comp",1:H)))
  ## weight matrix
  PP = matrix(0, p, H, dimnames=list(colnames(Xs), paste0("Comp",1:H)))
  ## R2X
  SSX = sum(Xs^2)
  R2X = matrix(0, 2, H, dimnames=list(c("pct","cum"), paste0("Comp",1:H)))
  ### SVD
  svd.X = svd(Xs, nu=H, nv=H)
  UU = svd.X$u
  VV = svd.X$v

  h = 1; Xh = Xs
  while(h<=H)
  {
    v.old = VV[, h, drop=F]
    u.old = UU[, h, drop=F]
    u.diff = u.old
    iter=1
    
    uv.h = .sparse.pca(u=u.old, v=v.old, Xh=Xh, lamda=lamda, tol=tol, max.iter=max.iter)
    u.new = uv.h$u.new
    v.new = uv.h$v.new
    
    Xhat = u.new %*% t(v.new)
    R2X["pct", h] = sum(Xhat^2) / SSX
    R2X["cum", h] = ifelse(h==1, R2X["pct", h], R2X["cum", h-1]+R2X["pct", h])
    # Deflation
    Xh = Xh - Xhat
    TT[,h] = u.new; PP[,h] = v.new
    h = h + 1
  }
 
  res.spca=structure(list(TT=TT, PP=PP, R2X=R2X, keep.x=keep.x, scale.x=scale.x, H=H), class="spca")
  return(invisible(res.spca))
}


.soft.threshold=function(input, lamda){
    # input=as.vector(input)
    if(lamda==0) return(input)
    cutoff=quantile(abs(input), lamda)
    output=abs(input)-cutoff
    output[output<=.Machine$double.eps]=0
    return(sign(input)*output)
}  
  
.sparse.pca=function(u, v, Xh, lamda, tol, max.iter){
  v.old = v
  u.old = u
  
  u.diff = u.old
  iter = 1
  
  p = nrow(v.old); n = nrow(u.old)
  repeat{
    #### Lasso: soft threshold
    v.new = matrix(.soft.threshold(t(Xh) %*% u.old, lamda), p, 1, dimnames=list(colnames(Xh), NULL))
    u.new = Xh %*% v.new
    u.new = u.new / sqrt(drop(t(u.new) %*% u.new))
    u.diff = u.old - u.new
    
    if(sum(u.diff^2) < tol | iter > max.iter) return(list(u.new = u.new, v.new = v.new))
    
    v.old = v.new
    u.old = u.new
    iter = iter + 1
  } 
} 
  
#TT=TT, PP=PP, R2X=R2X, keep.x=keep.x, scale.x=scale.x, H=H


vips.spca=function(spcaObj){
  PP=spcaObj$PP; p=nrow(PP); h=ncol(PP)
  R2X=spcaObj$R2[,1:h, drop=F]
  total.R2X=spcaObj$R2["cum",h]
  P2=PP^2 ### sum(P2)==1
  vips=sqrt(apply(P2%*%R2X, 1, sum)/total.R2X)*sqrt(p)
  return(vips)
}

