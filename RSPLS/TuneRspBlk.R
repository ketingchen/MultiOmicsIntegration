
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
       "/work/LAS/myn-lab/KChen/RSPLS/DataToBlocks.R",
       "/work/LAS/myn-lab/KChen/BlockPC.RF/sPCA.R")
for(pa in path) source(pa)

cat("Load input arguments","\n")
args = commandArgs(trailingOnly=T)
lab=c("-x","-y","-c","-slv","-elv","-m","-spca","-epca","-snF","-enF","-spls","-epls","-r","-k","-o","-nt")
idx=sapply(lab, function(ll)ifelse(sum(grepl(ll,args))==1,which(grepl(ll,args)),NA))
args1=rep("",length(lab));names(args1)=lab
for(i in 1:length(idx)){
        if(!is.na(idx[i])){
        args1[i]=gsub(names(idx)[i],"",args[idx[i]])} else args1[i]=NA
}
args=args1
print(args)

cat("Load input files:")
cat("x_")
x=read.csv(args["-x"],row.names=1)
cat("y_")
y=read.csv(args["-y"],row.names=1)

pc=args["-c"]

if(!is.na(as.numeric(pc))) y=y[,as.numeric(pc),drop=F]
if(pc=="rel") y=y[,c("rel.hc", "rel.fa")]
if(pc=="conc") y=y[,c("total","hc","fa")]
print(head(y))

if(is.na(args["-m"])) cls=rep(0, ncol(x)) else{
        synt.cls=read.csv(args["-m"], row.names=1)
        cls=synt.cls[,"Gene_Cluster"]
        cls=gsub("Cluster ","",cls)
        cls[cls=="Unassigned"]=0
        cls=as.numeric(cls)
        names(cls)=rownames(synt.cls)
}


if(is.na(args["-elv"])){
 if(is.na(args["-slv"])) lv=1 else lv=seq(as.numeric(args["-slv"]),ncol(y), 1)
} else{
 if(is.na(args["-slv"])) args["-slv"]=1
 lv=seq(as.numeric(args["-slv"]), as.numeric(args["-elv"]), 1)
} 

if(is.na(args["-spca"])|is.na(args["-epca"])) lamda.pca=0 else lamda.pca=seq(as.numeric(args["-spca"])/10,as.numeric(args["-epca"])/10,0.2)
if(is.na(args["-spls"])|is.na(args["-epls"])) lamda.rspls=0 else lamda.rspls=seq(as.numeric(args["-spls"])/10,as.numeric(args["-epls"])/10,0.2)
if(is.na(args["-snF"])|is.na(args["-enF"])) nFeatures=3 else nFeatures=seq(as.numeric(args["-snF"]),as.numeric(args["-enF"]),1)
 
#nFeatures=round(ncol(x)/10)


cat(lv,"\n")
#cat(nFeatures,"\n")
cat(lamda.pca,"\n")
cat(lamda.rspls,"\n")

cat("Tuning...","\n")
nt=as.numeric(args["-nt"])
tune=tuneBlk.cv.rspls(X=x, Y=y, cluster=cls, lv=lv, lamda.pca=lamda.pca, lamda.rspls=lamda.rspls, nFeatures=nFeatures, cv.fold=as.numeric(args["-k"]), cv.repeats=as.numeric(args["-r"]),nThreads=nt, plot=F)

if(!is.na(args["-o"])){
	outFile=paste0(args["-o"],".txt")
	write.table(tune,outFile)
}else print(tune)

cat("Job Finished", "\n")
