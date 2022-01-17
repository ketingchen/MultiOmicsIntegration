
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
require(randomForest)
require(doParallel)

path=c("/work/LAS/myn-lab/KChen/BlockPC.RF/updateX.R",
       "/work/LAS/myn-lab/KChen/BlockPC.RF/sPCA.R",
       "/work/LAS/myn-lab/KChen/BlockPC.RF/blockpc.RF.R",
       "/work/LAS/myn-lab/KChen/BlockPC.RF/variable.importance.R")
for(pa in path) source(pa)

args = commandArgs(trailingOnly=T)
lab=c("-x","-y","-b","-m","-c","-o","-r","-nt")
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
pc=as.numeric(args["-c"])

y=y[,pc]

synt.cls=read.csv(args["-m"], row.names=1)
cls=as.numeric(synt.cls[,"Gene_Cluster"])
names(cls)=rownames(synt.cls)
if(sum(cls==0)==0) cls[cls==max(cls)]=0

bprf=block.pc.rf(X=x, Y=y, lv=1, cluster=cls,quality.check=F, package="randomForest")

nt=as.numeric(args["-nt"])

cat("Calculate variable importance...","\n")
gene.rank=importance(bprf, idx.x=1:ncol(bprf$X.block),permutation=pp, nThreads=nt)
gene.rank=gene.rank[order(gene.rank[,"dR2"], gene.rank[,"dRMSE"], gene.rank[,"dAccuracy"], decreasing=T),]

r2=bprf$r2.forest
rmse=bprf$rmse.forest
accuracy=mean(bprf$accuracy)

if(!is.na(args["-o"])){
        outFile=paste0(args["-o"],".txt")
        write.table(gene.rank,outFile)
}
if(!is.na(args["-r"]))
{
        reportFile=paste0(args["-r"],".report")
        cat("BPRF:R2=",r2,"\n","RMSE=",rmse,"\n","Accuracy=",accuracy,"\n",file=reportFile)
}


cat("Job Finished", "\n")
