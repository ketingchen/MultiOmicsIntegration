## Instructions for querying the association between metabolome compositions and gene co-expression clusters using a combined approach of WGCNA and random forest

### 1. Input files needed:
1. a single column file that represent the composition of the entire metabolome using a single column of values, e.g., tSNE, first principal component, etc.<br>
2. Gene expression profiles, with each row representing a sample id, and each column represent represent a transcript. Please note that the sample info is not included in this file. <br>

### 2. WGCNA
This instruction only focuses on generating the co-expression clusters and querying the association between co-expressed clusters and the summarized metabolome composition scores (input file 1). Please refer to R package "Rtsne" or R functions prcomp() and princomp() for generating the single-value metabolome composition file (i.e., input file 1).<br>

2.1. Before start, selection the optimal thresholding power value to adjacency matrix. Please visit https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#46_Determine_parameters_for_WGCNA for details.<br>
 <br> 
2.2. Data cleaning. Please refer to function goodGenes() in R/WGCNA package, or nearZeroVar() in R/caret package to remove transcripts with too many missing entries or transcripts only expressed in a very small proportion of the samples (e.g., <5%).<br>
 <br>
2.3. WGCNA to generate co-expressed gene clusters. Please use R script "wgcna_methodA.R". The dataset used in thisscript (i.e., sample dataset) is available at /lss/research/myn-lab/Rupam/WGCNA/fpkm_clean_for_wgcna.csv <br>
 <br>
2.4. The output file to be used in the next step is "wgcna_cluster_membership_method_A.csv". This file is generated based on the sample dataset, and is available at /lss/research/myn-lab/Rupam/WGCNA/wgcna_cluster_membership_method_A.csv <br>

### 3. Random forest
Random forest regression is next constructed to query the association between co-expression gene clusters and single-value metabolome compositions.<br>

#### steps 3.1. and 3.2 are optional. Please use these two steps if you want to incorporate a lasso penality within each co-expressed clusters
##### 3.1. (optional) Generate the eigenvalue for each co-expressed gene cluster. Use DataToBlocks.R as demonstrated below (the sample dataset is available at: /lss/research/myn-lab/Rupam/BlockPC.RF/fpkm_for_cls.csv
      clusters <- read.csv("wgcna_cluster_membership_method_A.csv", row.names=1)
      cls <- as.numeric(cls[,"Gene_Cluster"])
      names(cls) <- rownames(clusters)
      if(sum(cls==0)==0) cls[cls==max(cls)] <- 0
      source("DataToBlocks.R")
      fpkm <- read.csv("fpkm_for_cls.csv", header=T, stringAsFactors=F, row.names=1)
      fpkm_new <- DataToBlocks(X=fpkm, cluster=cls)$Xnew

##### 3.2. (optional) Please use this step to find the optimal tuning parameters for lasso penalty. use tuneBPRF.R as demonstrated below:
     #!/bin/bash
     #SBATCH --nodes=1 
     #SBATCH --cpus-per-task=32 
     #SBATCH --mem=256G 
     #SBATCH --time=2-02:30:02 
     #SBATCH --partition=biocrunch
     #SBATCH --job-name="tuneECW" 
     #SBATCH --constraint=AVX2

     module load r/3.6.0-py2-fupx2uq
     module load r-randomforest/4.6-12-py2-r3.4-2biqljo
     module load r-doparallel/1.0.11-py2-r3.4-tzdczfv

     cd /work/LAS/myn-lab/Rupam/BlockPC.RF
     Rscript /work/LAS/myn-lab/KChen/BlockPC.RF/TuneBPRF.R -xfpkm_for_cls.csv -y../resp_tsne_rseq.csv -mwgcna_cluster_membership_method_A.csv -c2 -k4 -r3 -otune_Blk_tsne_ecw -nt32
                     
##### 3.3. Query the association between metabolome compositions and co-expressed gene clusters. Use RunBPRF.R as demonstrated below (the sample dataset for the response variables is available at: /lss/research/myn-lab/Rupam/resp_tsne_rseq.csv). The example script does NOT include lasso penalty within each co-expressed cluster.
     #!/bin/bash
     #SBATCH --nodes=1 
     #SBATCH --cpus-per-task=32 
     #SBATCH --mem=256G 
     #SBATCH --time=2-02:30:02 
     #SBATCH --partition=biocrunch
     #SBATCH --job-name="tuneECW" 
     #SBATCH --constraint=AVX2

     module load r/3.6.0-py2-fupx2uq
     module load r-randomforest/4.6-12-py2-r3.4-2biqljo
     module load r-doparallel/1.0.11-py2-r3.4-tzdczfv

     cd /work/LAS/myn-lab/Rupam/BlockPC.RF
     Rscript /work/LAS/myn-lab/KChen/BlockPC.RF/TuneBPRF.R -xfpkm_for_cls.csv -y../resp_tsne_rseq.csv -mwgcna_cluster_membership_method_A.csv -c2 -k4 -r3 -otune_Blk_tsne_ecw -nt32

    
      

