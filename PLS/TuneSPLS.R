
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
path=c("/work/LAS/myn-lab/KChen/SPLS/sPLS.R")
for(pa in path) source(pa)

cat("Load input arguments","\n")
args = commandArgs(trailingOnly=T)
lab=c("-x","-y","-c","-slv","-elv","-spls","-epls","-r","-k","-o")
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
cat("y_\n")
y=read.csv(args["-y"],row.names=1)

pc=args["-c"]

if(!is.na(as.numeric(pc))) y=y[,as.numeric(pc),drop=F]
if(pc=="rel") y=y[,c("rel.hc", "rel.fa")]
if(pc=="conc") y=y[,c("total","hc","fa")]
print(head(y))

if(is.na(args["-elv"])){
 if(is.na(args["-slv"])) lv=1 else lv=seq(as.numeric(args["-slv"]),ncol(y), 1)
} else{
 if(is.na(args["-slv"])) args["-slv"]=1
 lv=seq(as.numeric(args["-slv"]), as.numeric(args["-elv"]), 1)
}

if(is.na(args["-spls"])|is.na(args["-epls"])) lamda1=0 else lamda1=seq(as.numeric(args["-spls"])/10,as.numeric(args["-epls"])/10,0.1)

cat("Tuning...","\n")
tune=tune.spls(X=x, Y=y, lv=lv, lamda1=lamda1, cv.fold=as.numeric(args["-k"]), cv.repeats=as.numeric(args["-r"]), plot=F)

if(!is.na(args["-o"])){
	outFile=paste0(args["-o"],".txt")
	write.table(tune,outFile)
}else print(tune)

cat("Job Finished", "\n")
