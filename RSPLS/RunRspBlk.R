
# Sample:
# RunSPLS.R(0) -xpred_rpkm_14.rds -yresp_pc_14.rds -b10000 -c1 -oReturn -rY -nt16
#     -x: predictors
#     -y: responses
#     -b: permutation runs
#     -c: the index of column to be used as the  response variable
#     -nt: nThreads
#     -o: path to write the output (optional)
#     -r: path to write the model evaluation report (optional)
#   
require(doSNOW)

path=c("/work/LAS/myn-lab/KChen/RSPLS/sPLS.R",
       "/work/LAS/myn-lab/KChen/RSPLS/randomSPLS.multi.R",
       "/work/LAS/myn-lab/KChen/BlockPC.RF/sPCA.R",
       "/work/LAS/myn-lab/KChen/RSPLS/DataToBlocks.R")
for(pa in path) source(pa)

args = commandArgs(trailingOnly=T)
lab=c("-x","-y","-b","-m", "-lv","-c","-lpca","-lpls","-o","-r","-nt")
idx=sapply(lab, function(ll)ifelse(sum(grepl(ll,args))==1,which(grepl(ll,args)),NA))
args1=rep("",length(lab));names(args1)=lab
for(i in 1:length(idx)){
        if(!is.na(idx[i])){
        args1[i]=gsub(names(idx)[i],"",args[idx[i]])} else args1[i]=NA
}
args=args1

x=read.csv(args["-x"],row.names=1)
y=read.csv(args["-y"],row.names=1)

pp=as.numeric(args["-b"])
pc=args["-c"]

if(!is.na(as.numeric(pc))) y=y[,as.numeric(pc),drop=F]
if(pc=="rel") y=y[,c("rel.hc", "rel.fa")]
if(pc=="conc") y=y[,c("total","hc","fa")]

if(is.na(args["-m"])) cls=rep(0, ncol(x)) else{
        synt.cls=read.csv(args["-m"], row.names=1)
        cls=synt.cls[,"Gene_Cluster"]
        cls=gsub("Cluster ","",cls)
        cls[cls=="Unassigned"]=0
        cls=as.numeric(cls)
        names(cls)=rownames(synt.cls)
}

lv=as.numeric(ifelse(is.na(args["-lv"]), 1, as.numeric(args["-lv"])))
lamda.pca=as.numeric(ifelse(is.na(args["-lpca"]), 0, as.numeric(args["-lpca"])))/10
lamda.rspls=as.numeric(ifelse(is.na(args["-lpls"]), 0, as.numeric(args["-lpls"])))/10

xnew=DataToBlocks(x, cluster=cls, lamda=rep(lamda.pca, max(cls)))$Xnew

nt=as.numeric(args["-nt"])
print(dim(xnew))
print(dim(y))

rsp=randomSPLS(X=xnew, Y=y, lv=lv, lamda1=lamda.rspls,quality.check=F, verbose=T, nThreads=nt)

cat("Calculate variable importance...","\n")

pvalue=permute(rsp, PP=pp, nThreads=nt)
    
padj=pvalue[,2]
padj[pvalue[,1]>1]=p.adjust(padj[pvalue[,1]>1],method="BH")
    
r2=rsp$r2.forest
rmse=rsp$rmse.forest
accuracy=rsp$accuracy.forest

gene.rank=cbind(importance=rsp$importance,pvalue, padj=padj)
gene.rank=gene.rank[order(gene.rank[,"importance"],decreasing=T),]

if(!is.na(args["-o"])){
	outFile=paste0(args["-o"],".txt")
	write.table(gene.rank,outFile)
}
if(!is.na(args["-r"]))
{
        reportFile=paste0(args["-r"],".report")
        cat("RsplsBlk:R2=",r2,"\n","RMSE=",rmse,"\n","Accuracy=",accuracy,"\n",file=reportFile)         
}  

cat("Job Finished", "\n")
