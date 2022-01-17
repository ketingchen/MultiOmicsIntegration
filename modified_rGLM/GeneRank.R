gList=function(RDS,B=10000){
  count=0
  gl=function(RR){
   size=length(RR)
	for(i in 1:size){
	rds=readRDS(RR[[i]])
	count=count+rds[,"count"]
	if(i==1)g0=rds[,"g0"]
	rm(rds)
	}
  pval=(1+count)/(1+B)
  names(pval)=names(count)=names(g0)
  return(cbind(g0, count, pval))
  }
 lapply(RDS, function(RR)gl(RR))
}

compList=function(totalGenes,geneList){
   size=length(geneList)
   m=matrix(0,nrow=length(totalGenes),ncol=length(c("g0","pval"))*size,
     dimnames=list(totalGenes,paste0(c("g0_PC","pval_PC"),rep(1:size,each=2))))
   for(i in 1:size){
	gene=geneList[[i]]
	g0=paste0("g0_PC",i);pval=paste0("pval_PC",i)
	m[rownames(gene),g0]=gene[,"g0"]
        m[rownames(gene),pval]=gene[,"pval"]
        m[m[,pval]==0,pval]=1
   }
  return(m)
}

gRank.default=function(m,geneListSize=2){
  size=geneListSize
  pvals=paste0("pval_PC",1:size)
  PVAL=apply(m[,pvals],1, min)
  g0s=paste0("g0_PC",1:size)
  G0=apply(m[,g0s],1,max)
  PADJ=PVAL
  PADJ[G0>0]=p.adjust(PADJ[G0>0],method="BH")
  RANK=rep(0,length(PVAL));names(RANK)=names(PVAL)
  RANK[PVAL<0.05]=1
  if(sum(PADJ<0.05)==0)	RANK[PVAL==min(PVAL)]=2
  else RANK[PADJ<0.05]=2
  
  gr=cbind(G0,PVAL,PADJ,RANK)
  gr=gr[order(gr[,"RANK"],gr[,"G0"],decreasing=T),]
  return(gr)
}

 gRank=function(totalGenes, RDS, B=10000){
  geneList=gList(RDS,B)
  geneMatrix=compList(totalGenes, geneList)
  return(gRank.default(geneMatrix,length(geneList)))	
}
