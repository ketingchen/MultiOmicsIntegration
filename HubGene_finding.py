# This script is created according to Das et al.(2017) DOI:10.1371/journal.pone.0169605
import pandas as pd
import sys
sys.path.append("/Users/kchen/Box/MAIZE SEEDLING MANUSCRIPT_Rupam and Keting/Keting/HubGenes")
from FindHub import main, edgelist, findHub

cluster1 = pd.read_csv("Data/fpkm_cluster1_brown.csv", index_col=0)
cluster19 = pd.read_csv("Data/fpkm_cluster19_white.csv", index_col=0)
cluster2 = pd.read_csv("Data/fpkm_cluster2_royalblue.csv", index_col=0)
cluster6 = pd.read_csv("Data/fpkm_cluster6_black.csv", index_col=0)
cluster8 = pd.read_csv("Data/fpkm_cluster8_darkorange.csv", index_col=0)
cluster9 = pd.read_csv("Data/fpkm_cluster9_blue.csv", index_col=0)
cluster10 = pd.read_csv("Data/fpkm_cluster10_yellow.csv", index_col=0)
cluster11 = pd.read_csv("Data/fpkm_cluster11_orange.csv", index_col=0)
cluster12 = pd.read_csv("Data/fpkm_cluster12_paleturqoise.csv", index_col=0)
cluster14 = pd.read_csv("Data/fpkm_cluster14_salmon.csv", index_col=0)
cluster18 = pd.read_csv("Data/fpkm_cluster18_orangered4.csv", index_col=0)


modules = pd.read_csv("Data/wgcna_ME_method_A.csv", index_col=0)
traits = pd.read_csv("Data/resp_tsne_rseq.txt", sep=" ")



p1 = main(cluster1, beta=1, S=1000)
p19 = main(cluster19, beta=1, S=1000)
p2 = main(cluster2, beta=1, S=1000)
p6 = main(cluster6, beta=1, S=1000)
p8 = main(cluster8, beta=1, S=1000)
p9 = main(cluster9, beta=1, S=1000)
p10 = main(cluster10, beta=1, S=1000)
p11 = main(cluster11, beta=1, S=1000)
p12 = main(cluster12, beta=1, S=1000)
p14 = main(cluster14, beta=1, S=1000)
p18 = main(cluster18, beta=1, S=1000)

p1.to_csv("Cluster1/Cluster1_pvalues_connection_degree.csv")
p2.to_csv("OtherClusters/Cluster2_pvalues_connection_degree.csv")
p6.to_csv("OtherClusters/Cluster6_pvalues_connection_degree.csv")
p8.to_csv("OtherClusters/Cluster8_pvalues_connection_degree.csv")
p9.to_csv("OtherClusters/Cluster9_pvalues_connection_degree.csv")
p10.to_csv("OtherClusters/Cluster10_pvalues_connection_degree.csv")
p11.to_csv("OtherClusters/Cluster11_pvalues_connection_degree.csv")
p12.to_csv("OtherClusters/Cluster12_pvalues_connection_degree.csv")
p14.to_csv("OtherClusters/Cluster14_pvalues_connection_degree.csv")
p18.to_csv("OtherClusters/Cluster18_pvalues_connection_degree.csv")
p19.to_csv("Cluster19/Cluster19_pvalues_connection_degree.csv")


##### Reduce the giant correlation network

p1 = pd.read_csv("Cluster1/Cluster1_pvalues_connection_degree.csv", index_col=0)
p2 = pd.read_csv("OtherClusters/Cluster2_pvalues_connection_degree.csv", index_col=0)
p6 = pd.read_csv("OtherClusters/Cluster6_pvalues_connection_degree.csv", index_col=0)
p8 = pd.read_csv("OtherClusters/Cluster8_pvalues_connection_degree.csv", index_col=0)
p9 = pd.read_csv("OtherClusters/Cluster9_pvalues_connection_degree.csv", index_col=0)
p10 = pd.read_csv("OtherClusters/Cluster10_pvalues_connection_degree.csv", index_col=0)
p11 = pd.read_csv("OtherClusters/Cluster11_pvalues_connection_degree.csv", index_col=0)
p12 = pd.read_csv("OtherClusters/Cluster12_pvalues_connection_degree.csv", index_col=0)
p14 = pd.read_csv("OtherClusters/Cluster14_pvalues_connection_degree.csv", index_col=0)
p18 = pd.read_csv("OtherClusters/Cluster18_pvalues_connection_degree.csv", index_col=0)
p19 = pd.read_csv("Cluster19/Cluster19_pvalues_connection_degree.csv", index_col=0)



edgelist1 = edgelist(cluster1, beta=1)
edgelist19 = edgelist(cluster19, beta=1)
edgelist1.to_csv('Cluster1_EdgeList.csv')
edgelist19.to_csv('Cluster19_EdgeList.csv')

h1 = findHub(p1, cluster1, traits, modules['MEbrown'], corr_cut=0.7)[3]
h19 = findHub(p19, cluster19, traits, modules['MEwhite'], corr_cut=0.75)[3]
k1 = cluster1.columns.to_list()
hb1 = [k in h1 for k in k1]
k19 = cluster19.columns.to_list()
hb19 = [k in h19 for k in k19]

node1 = pd.DataFrame({'Genes':k1, 'Hub':hb1})
node19 = pd.DataFrame({'Genes':k19, 'Hub':hb19})
node1.to_csv('Cluster1_NodeTable.csv')
node19.to_csv('Cluster19_NodeTable.csv')

reduce_idx1 = list(set([i for i, e in enumerate(edgelist1['Gene_A']) if e in h1] + [j for j, e in enumerate(edgelist1['Gene_B']) if e in h1]))
edgelist1_reduce = edgelist1.iloc[reduce_idx1, :]
#### Add another layer to screen for correlations with abs value >0.85
#mask = (abs(edgelist1_reduce['corr'])>0.8)
#edgelist1_reduce = edgelist1_reduce[mask]
node1_reduce = list(set([a for a in edgelist1_reduce['Gene_A']] + [b for b in edgelist1_reduce['Gene_B']]))
hbr1 = [k in h1 for k in node1_reduce]
edghbr1 = [a in h1 and b in h1 for a, b in zip(edgelist1_reduce['Gene_A'], edgelist1_reduce['Gene_B'])]
edgelist1_reduce['edge_hub'] = edghbr1
node1_reduce = pd.DataFrame({'Nodes':node1_reduce, 'Hub':hbr1})
edgelist1_reduce.to_csv("Cluster1_EdgeList_Reduce_hub15.csv")
node1_reduce.to_csv("Cluster1_NodeTable_Reduce_hub15.csv")

reduce_idx19 = list(set([i for i, e in enumerate(edgelist19['Gene_A']) if e in h19] + [j for j, e in enumerate(edgelist19['Gene_B']) if e in h19]))
edgelist19_reduce = edgelist19.iloc[reduce_idx19, :]
#### Add another layer to screen for correlations with abs value >0.85
mask = (abs(edgelist19_reduce['corr'])>0.7)
edgelist19_reduce = edgelist19_reduce[mask]
node19_reduce = list(set([a for a in edgelist19_reduce['Gene_A']] + [b for b in edgelist19_reduce['Gene_B']]))
hbr19 = [k in h19 for k in node19_reduce]
node19_reduce = pd.DataFrame({'Nodes':node19_reduce, 'Hub':hbr19})
edgelist19_reduce.to_csv("Cluster19_EdgeList_Reduce_hub23_edgeCutPct70.csv")
node19_reduce.to_csv("Cluster19_NodeTable_Reduce_hub23.csv")

##### The connection within the hubs
nethub1 = edgelist1[edgelist1['Gene_A'].isin(h1) & edgelist1['Gene_B'].isin(h1)]
nethub19 = edgelist19[edgelist19['Gene_A'].isin(h19) & edgelist19['Gene_B'].isin(h19)]
nethub1.to_csv("Cluster1_EdgeList_Hub140_only.csv")
nethub19.to_csv("Cluster19_EdgeList_Hub185_only.csv")

deg1_GeneA = dict([[h, edgelist1['Gene_A'].tolist().count(h)] for h in h1])
deg1_GeneB = dict([[h, edgelist1['Gene_B'].tolist().count(h)] for h in h1])
deg1 = [deg1_GeneA[h] + deg1_GeneB[h] for h in h1]
node1_hub = pd.DataFrame({'Genes':h1, 'Degrees':deg1})
node1_hub.to_csv("Cluster1_NodeTable_Hub140_only.csv")

deg19_GeneA = dict([[h, edgelist19['Gene_A'].tolist().count(h)] for h in h19])
deg19_GeneB = dict([[h, edgelist19['Gene_B'].tolist().count(h)] for h in h19])
deg19 = [deg19_GeneA[h] + deg19_GeneB[h] for h in h19]
node19_hub = pd.DataFrame({'Genes':h19, 'Degrees':deg19})
node19_hub.to_csv("Cluster19_NodeTable_Hub185_only.csv")

#### Find the hubs of other clusters
h1 = findHub(p1, cluster1, traits, modules['MEbrown'], corr_cut=0.7)[3]
h2 = findHub(p2, cluster2, traits, modules['MEroyalblue'], corr_cut=0.7)[3]
h6 = findHub(p6, cluster6, traits, modules['MEblack'], corr_cut=0.6, whichTrait="cutin")[3]
h8 = findHub(p8, cluster8, traits, modules['MEdarkorange'], corr_cut=0.5)[3]
h9 = findHub(p9, cluster9, traits, modules['MEblue'], corr_cut=0.65)[3]
h10 = findHub(p10, cluster10, traits, modules['MEyellow'], corr_cut=0.55, whichTrait="cutin")[3]
h11 = findHub(p11, cluster11, traits, modules['MEorange'], corr_cut=0.45, whichTrait="cutin")[3]
h12 = findHub(p12, cluster12, traits, modules['MEpaleturquoise'], corr_cut=0.3)[3]
h14 = findHub(p14, cluster14, traits, modules['MEsalmon'], corr_cut=0.3)[3]
h18 = findHub(p18, cluster18, traits, modules['MEorangered4'], corr_cut=0.65)[3]
h19 = findHub(p19, cluster19, traits, modules['MEwhite'], corr_cut=0.75)[3]

candidates = pd.read_csv("../Candidates/candidates_total.txt", sep=' ')['Zm_ID'].to_list()
hub_file1 = pd.DataFrame({'Hub':h1, 'Candidates':[h in candidates for h in h1]})
hub_file2 = pd.DataFrame({'Hub':h2, 'Candidates':[h in candidates for h in h2]})
hub_file6 = pd.DataFrame({'Hub':h6, 'Candidates':[h in candidates for h in h6]})
hub_file8 = pd.DataFrame({'Hub':h8, 'Candidates':[h in candidates for h in h8]})
hub_file9 = pd.DataFrame({'Hub':h9, 'Candidates':[h in candidates for h in h9]})
hub_file10 = pd.DataFrame({'Hub':h10, 'Candidates':[h in candidates for h in h10]})
hub_file11 = pd.DataFrame({'Hub':h11, 'Candidates':[h in candidates for h in h11]})
hub_file12 = pd.DataFrame({'Hub':h12, 'Candidates':[h in candidates for h in h12]})
hub_file14 = pd.DataFrame({'Hub':h14, 'Candidates':[h in candidates for h in h14]})
hub_file18 = pd.DataFrame({'Hub':h18, 'Candidates':[h in candidates for h in h18]})
hub_file19 = pd.DataFrame({'Hub':h19, 'Candidates':[h in candidates for h in h19]})

hub_file1.to_csv("HubGenes_Cluster1.txt", sep="\t", index=False)
hub_file2.to_csv("HubGenes_Cluster2.txt", sep="\t", index=False)
hub_file6.to_csv("HubGenes_Cluster6.txt", sep="\t", index=False)
hub_file8.to_csv("HubGenes_Cluster8.txt", sep="\t", index=False)
hub_file9.to_csv("HubGenes_Cluster9.txt", sep="\t", index=False)
hub_file10.to_csv("HubGenes_Cluster10.txt", sep="\t", index=False)
hub_file11.to_csv("HubGenes_Cluster11.txt", sep="\t", index=False)
hub_file12.to_csv("HubGenes_Cluster12.txt", sep="\t", index=False)
hub_file14.to_csv("HubGenes_Cluster14.txt", sep="\t", index=False)
hub_file18.to_csv("HubGenes_Cluster18.txt", sep="\t", index=False)
hub_file19.to_csv("HubGenes_Cluster19.txt", sep="\t", index=False)