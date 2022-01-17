
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
lab=c("-x","-y","-c","-m","-k","-r","-o","-nt")
idx=sapply(lab, function(ll)ifelse(sum(grepl(ll,args))==1,which(grepl(ll,args)),NA))
args1=rep("",length(lab));names(args1)=lab
for(i in 1:length(idx)){
	if(!is.na(idx[i])){
	args1[i]=gsub(names(idx)[i],"",args[idx[i]])} else args1[i]=NA
}
args=args1

x=read.csv(args["-x"],row.names=1)
y=read.csv(args["-y"],row.names=1)

pc=as.numeric(args["-c"])

y=y[,pc]

synt.cls=read.csv(args["-m"], row.names=1)
cls=as.numeric(synt.cls[,"Gene_Cluster"]) 
names(cls)=rownames(synt.cls)
if(sum(cls==0)==0) cls[cls==max(cls)]=0

cv.fold=as.numeric(args["-k"])
cv.repeats=as.numeric(args["-r"])

tune=tuneLamda(X=x, Y=y, cluster=cls, lamda=seq(0, 0.8, 0.2), plot=F, cv.fold=cv.fold, cv.repeats=cv.repeats)

write.table(tune, args["-o"])
cat("Job Finished", "\n")
