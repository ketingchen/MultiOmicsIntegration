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
