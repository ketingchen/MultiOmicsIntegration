
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

path=c("/work/LAS/myn-lab/KChen/SPLS/sPLS.R")
for(pa in path) source(pa)

args = commandArgs(trailingOnly=T)
lab=c("-x","-y","-b","-c","-lv","-lx","-ly","-o","-nt")
idx=sapply(lab, function(ll)ifelse(sum(grepl(ll,args))==1,which(grepl(ll,args)),NA))
args1=rep("",length(lab));names(args1)=lab
for(i in 1:length(idx)){
        if(!is.na(idx[i])){
        args1[i]=gsub(names(idx)[i],"",args[idx[i]])} else args1[i]=NA
}
args=args1

print(args)
x=read.csv(args["-x"],row.names=1)
y=read.csv(args["-y"],row.names=1)

pp=as.numeric(args["-b"])
pc=args["-c"]

if(!is.na(as.numeric(pc))) y=y[,as.numeric(pc),drop=F]
if(pc=="rel") y=y[,c("rel.hc", "rel.fa")]
if(pc=="conc") y=y[,c("total","hc","fa")]
print(head(y))

lv=ifelse(is.na(args["-lv"]), 1, args["-lv"])
lv=as.numeric(lv)
lamda1=as.numeric(ifelse(is.na(args["-lx"]), 0, args["-lx"]))/10
lamda2=as.numeric(ifelse(is.na(args["-ly"]), 0, args["-ly"]))/10

nt=as.numeric(args["-nt"])

print(dim(x))
print(dim(y))
print(lv)
print(lamda1)
print(lamda2)

# sPLS=function(X, Y, lv=1, max.iter=500, tol=1e-06, lamda1=0, lamda2=0, scale=T, verbose=T, quality.check=T)
sp=sPLS(X=x, Y=y, lv=lv, lamda1=lamda1, lamda2=lamda2, scale=T, quality.check=F)

cat("Calculate p-values for variable importance","\n")
gene.rank=permute(sp, B=pp, nThreads=nt)
gene.rank=gene.rank[order(gene.rank[,"importance"],decreasing=T),]

cat("Permutation finished\n")
if(!is.na(args["-o"])){
	outFile=paste0(args["-o"],".txt")
	print(outFile)
	write.table(gene.rank,outFile)
}

cat("Job Finished", "\n")
