#### Instructions for querying the association between metabolome compositions and gene co-expression clusters using a combined approach of WGCNA and random forest

##### Input files needed:
1. a single column file that represent the composition of the entire metabolome using a single column of values, e.g., tSNE, first principal component, etc.<br>
2. Gene expression profiles, with each row representing a sample id, and each column represent represent a transcript <br>

##### WGCNA
This instruction only focuses on generating the co-expression clusters and querying the association between co-expressed clusters and the summarized metabolome composition scores (input file 1). Please refer to R package "Rtsne" or R functions prcomp() and princomp() for generating the single-value metabolome composition file (i.e., input file 1).<br>

1. Before start, selection the optimal thresholding power value to adjacency matrix. Please visit https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#46_Determine_parameters_for_WGCNA for details.
