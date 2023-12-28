####### RUPAM WGCNA ######

###########FOLLOWING WGCNA TUTORIAL STEP 1 - LOADING DATA AND REMOVING NOISE/OUTLIERS#######
require(WGCNA)

#### CHECK IF DATASET NEEDS TO BE TRANSCPOSED ####
##### READ.CSV FUNC - PROVIDE ABSOLUTE PATH NAME, HEADER = T, STRINGSASFACTORS = FALSE, AND ROW.NAMES IF DESIRED#####
#transpose so that genes are columns and samples are rows
fpkm = read.csv("/work/LAS/myn-lab/Rupam/WGCNA/fpkm_clean_for_wgcna.csv", header = T, stringsAsFactors = FALSE, row.names = 1)

#### CHECK WHERE NUMERIC DATA STARTS FOR DATASET  #######
#### cleaned_fpkm_seedling_data starts at col 4
#### rupam_fpkm_vips.csv starts at col 5 ####
exprData0 = fpkm
cat("Data loaded and expression data assigned to variable\n")

colnames = colnames(exprData0)
rownames = rownames(exprData0)

nGenes = ncol(exprData0)
nSamples = nrow(exprData0)
cat("preliminary data cleaning complete\n")

######## Enable parallel computing for WGCNA, with 32 CPUs ########
allowWGCNAThreads(nThreads=32)

###########FOLLOWING WGCNA TUTORIAL 2B - STEP BY STEP NETWORK CONSTRUCTION AND MODULE DETECTION#########
# Choose a set of soft-thresholding powers - soft thresholding power β to which co-expression similarity is raised to, i.e. (coexpression val)^(thresholding power)

#creating adjacency matrix 
#signed network analysis and with biweight midcorrelation 
#power 16 was chosen by examinging pickSoftThreshold function, powers 1-4  did not conform to normal Truncated R^2 values.
# 16 was first power with Truncated R^2 > 0.7
adjacency = adjacency(exprData0, type = "signed", corFnc = "bicor", power = 16)
cat("adjacency matrix created with soft power threshold\n")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
cat("TOM dissimilarity calculated\n")

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 30

# deepSplit has range of 0 to 4, and ths higher the value is, the more small cluster will be produced
dynamic_mods = cutreeDynamic(dendro=geneTree, distM=dissTOM,  method="hybrid", deepSplit = 2, minClusterSize = minModuleSize)
cat("gene cluster modules created\n")

#convert numeric labels into colors, and using that data merge modules that are above 0.15 (cutHeight = Maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging)
moduleColors = labels2colors(dynamic_mods)

merge = mergeCloseModules(exprData0, moduleColors, corFnc = "bicor", cutHeight = 0.15, verbose = 3)
merged_mods = data.frame(table(merge$colors))

#write.csv(merged_mods, "mods_gene_freq_methodB_pow14.csv")

mergedColors = merge$colors
mergedColors[mergedColors!="grey"] = paste0("Cluster_",mergedColors[mergedColors!="grey"])
mergedColors[mergedColors=="grey"] = "Unassigned"

mergedMEs = merge$newMEs
cat("merged ME color and eigengenes created")

#add columns that label eigengene sample genotype and tissue type
sf=read.table("/work/LAS/myn-lab/Rupam/FPKM_rupam_scr_sampleinfo.txt")
sf=sf[,1:2]

mergedMEs_labeled = data.frame(sf, mergedMEs[,1:ncol(mergedMEs)])

colors=colnames(mergedMEs_labeled)[-1:-2]
colors=gsub("ME","Cluster_",colors)
colors[colors=="Cluster_grey"]="Unassigned"

dict=data.frame(colors=colors, num=1:length(colors))
ModuleInNum=NULL
for(mm in mergedColors){
   nn=dict[dict$colors==mm, 2]
   ModuleInNum=c(ModuleInNum, nn)
}

write.csv(mergedMEs_labeled, "wgcna_ME_method_A.csv")
write.csv(dict, "Module_colorcode_idx_dictionary_methodA.csv")
#cat("Method B pow 14  merged eigengene list created")


####ONLY CALCULATE IF EIGENGENE DENDROGRAM IS REQUESTED ####
#mergedMEDiss = 1 - cor(mergedMEs)
#METree = hclust(as.dist(mergedMEDiss), method = "average")

#pdf("mergedME_dendrogram.pdf")
#plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#dev.off()
#cat("eigengene dendrogram created")

merged_module_summary = data.frame(Gene = colnames(exprData0), Module = mergedColors, ModuleInNum)

write.csv(merged_module_summary, "wgcna_cluster_membership_method_A.csv")
cat("summary table of genes created\n")


############SUMMARY FIGURES AND TABLES###########
#plot and save gene + color module dendrogram
#pdf("merged_gene_dendrogram_methodA.pdf")
#plotDendroAndColors(geneTree, cbind(moduleColors,mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
#       dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
#        guideHang = 0.05, main = "Gene dendrogram and merged module colors")
#dev.off()
#cat("gene dendrogram created")

