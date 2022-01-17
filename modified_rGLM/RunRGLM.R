# Rscript DimRdc_rglm.R(0) DataFile_pred(1) DataFile_resp(2) permRuns(3) PC(4)  outFile(5) wrtieReport(6) 

# Sample:
# DimRdc_rglm.R(0) -xpred_rpkm_14.rds -yresp_pc_14.rds -b1000 -c1 -oReturn.rds -rY
#     -x: predictors
#     -y: responses
#     -b: permutation runs
#     -c: the index of column to be used as the  response variable
#     -o: path to write the output (optional)
#     -r: path to write the model evaluation report (optional)

source("/work/LAS/myn-lab/KChen/RGLM/crossRGLM.R")
source("/work/LAS/myn-lab/KChen/PLS/PLS.R")

args = commandArgs(trailingOnly=T)
lab=c("-x","-y","-b","-c","-o","-r","-nt")

idx=sapply(lab, function(ll)ifelse(sum(grepl(ll,args))==1,which(grepl(ll,args)),NA))
args1=rep("",length(lab));names(args1)=lab
for(i in 1:length(idx)){
        if(!is.na(idx[i])){
        args1[i]=gsub(names(idx)[i],"",args[idx[i]])} else args1[i]=NA
}
args=args1

x=read.csv(args["-x"], row.names=1)
y=read.csv(args["-y"], row.names=1)

bb = as.numeric(args["-b"])
pc = as.numeric(args["-c"])

y=y[,pc]

plsm=pls.nipals(X=x, Y=y, lv=10, scale=T)
R2Y=sum(plsm$R2Y)
vips.pls=vips(plsm)
x.new=x[,vips.pls>1]

cores=as.numeric(args["-nt"])

cat("Number workers:",cores,"\n")

cat("Start permutation\n")

w = whole.rglm(X=x.new, Y=y, B=bb, numCores=cores)
  
cat("gene.rank: count=sapply(1:length(g0),function(j)sum(g1[j,]>=g0[j])); g0=genes0; pval=(1+count)/(1+bb)","\n")

g0=w$genes0
g1=w$genes1
estimate0=w$estimate0
rm(w)
cor.e=cor(estimate0)[1,2]
r2.e=summary(lm(estimate0))$r.squared

count=sapply(1:length(g0),function(j)sum(g1[j,]>=g0[j]))
names(count)=names(g0) 
pval=(1+count)/(1+bb)
padj=pval
padj[g0>=1]=p.adjust(padj[g0>=1],method="BH")

gene.rank=cbind(g0, count, pval, padj)

cat("Computation Finished","\n")

if(!is.na(args[5])){
        outFile=paste0(args["-o"],".rds")
	saveRDS(gene.rank, file=outFile)
}
if(!is.na(args[6]))
{
       reportFile=paste0(args["-r"],".report")
       cat("PLS:R2Y=",R2Y,"\n","Estimate corr=",cor.e,"\n","Estimate R2=",r2.e,"\n",file=reportFile)         
}  

cat("Job Finished", "\n")
