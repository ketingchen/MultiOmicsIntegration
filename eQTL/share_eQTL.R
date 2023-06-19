setwd("~/Box/SILKY/Keting/Projects/Omics_Integration/eQTL/Final_Figs_Tables_forManuscript")
### genomic link plot
library(circlize)
genetic_map <- read.csv("../Genetic_Map.csv")
b73v4 <- read.table("../B73_v4_gene_chr_position.bed", sep="\t")[358:44352,]
b73v4$V1 <- gsub("chr", "Chr ", b73v4$V1)
colnames(b73v4) <- c("chr", "start", "end", "gene")
b73v4$chr <- factor(b73v4$chr, levels=paste0("Chr ",1:10))

#### format for circular genome layout
chr <- paste0("Chr ", as.vector(genetic_map[,2]))
start <- end <- NULL
for(i in 1:nrow(genetic_map)){
	to_remove <- paste0(genetic_map[i,"Chr"], "_")
	start <- c(start, as.numeric(gsub(to_remove, "", genetic_map[i, "VariantID"])))
	end <- c(end, start[i]+1)
}

marker <- as.vector(genetic_map[,"cxxmxxxID"])
map <- data.frame(chr, start, end, marker)
map$chr <- factor(map$chr, levels=unique(chr))

snp <- data.frame(chr=map$chr, start=as.vector(map$start)-1, end=as.vector(map$end), value=as.vector(genetic_map[,"Kosambi_Distance"]))

#### combine b73v4 gene pos and snp pos to create a comprehensive genomic map
genome <- rbind(b73v4[,1:3], map[,1:3])


#### identify specific traits and eQTLs
candidate_eQTL <- read.csv("candidate_eQTL_genes_summary_c&d.csv")
candidate_eQTL <- candidate_eQTL[,c("qtl_chr", "qtl_start", "qtl_end", "peak_lod","gene")]

transcriptome_eQTL <- read.csv("transcriptome_eQTL_genes_summary_c&d.csv")
trait <- transcriptome_eQTL[,c(1:4, 8, 9)]
transcriptome_eQTL <- transcriptome_eQTL[,c("qtl_chr", "qtl_start", "qtl_end", "peak_lod","gene")]


cuticle_can <- unique(as.vector(read.csv("Regulatory_relationship.csv", header=F)$V1))
cuticle_can <- sort(c(cuticle_can, "Zm00001d023343", "Zm00001d027488", "Zm00001d027766", "Zm00001d036197", "Zm00001d043049", "Zm00001d044907", "Zm00001d048370", "Zm00001d051932"))
cuticle_can <- cuticle_can[-c(2, 14)]
cuticle_can_eQTL <- candidate_eQTL[candidate_eQTL$gene %in% cuticle_can,]

cuticle <- sort(unique(c(as.vector(read.table("../../Candidates/maize_lipid_orthologs")$V1), unique(as.vector(read.csv("Regulatory_relationship.csv", header=F)$V1)), "Zm00001d023343", "Zm00001d027488", "Zm00001d027766", "Zm00001d036197", "Zm00001d043049", "Zm00001d044907", "Zm00001d048370", "Zm00001d051932", "Zm00001d030343", "Zm00001d039094", "Zm00001d018455", "Zm00001d009608", "Zm00001d028241", "Zm00001d044579", "Zm00001d032948", "Zm00001d045660", "Zm00001d027904", "Zm00001d016438", "Zm00001d020206", "Zm00001d047931", "Zm00001d028406", "Zm00001d039053", "Zm00001d014168", "Zm00001d037328", "Zm00001d027766", "Zm00001d032728", "Zm00001d046444", "Zm00001d009354", "Zm00001d044580", "Zm00001d051787", "Zm00001d029350", "Zm00001d047190", "Zm00001d021200")))
cuticle_eQTL <- transcriptome_eQTL[transcriptome_eQTL$gene %in% cuticle,]


overlap <- function(start0, end0, start1, end1){
	### the 2nd eQTL locates in the upstream of the first eQTL
	if(end1 <= start0) return(0)
	### the 2nd eQTL locates in the downstream of the first eQTL
	if(start1 >= end0) return(0)
	### the rest scenarios are when the 2nd eQTL overlaps with the first, either in the upstream side or in the down stream side
	return(1)
}

share_eQTL <- function(query, query_chr, query_start, query_end, subject, annotation){
	share <- NULL
	for(i in 1:nrow(subject)){
		
		sub_chr <- subject[i,"qtl_chr"]
		sub_start <- subject[i,"qtl_start"]
		sub_end <- subject[i,"qtl_end"]
		sub <- subject[i,"gene"]
		
		if(query_chr == sub_chr){
			if(overlap(query_start, query_end, sub_start, sub_end)==1 & sub!=query){
			share <- rbind(share, data.frame(query, chr=query_chr, query_start, query_end, sub_start, sub_end, sub))
			}
		}
	}
	return(share)
}


#### finding candidates with overlapping eQTLs with the cuticle candidates
share_cuticle_can <- NULL
for(i in 1:nrow(cuticle_can_eQTL)){
	query <- cuticle_can_eQTL[i, "gene"]
	query_chr <- cuticle_can_eQTL[i, "qtl_chr"]
	query_start <- cuticle_can_eQTL[i, "qtl_start"]
	query_end <- cuticle_can_eQTL[i, "qtl_end"]
	share_cuticle_can <- rbind(share_cuticle_can, share_eQTL(query=query, query_chr=query_chr, query_start=query_start, query_end=query_end, subject=candidate_eQTL))
}
### remove genes with problometic eQTLs
to_remove <- c("Zm00001d018309", "Zm00001d048311")
share_cuticle_can <- share_cuticle_can[!share_cuticle_can$sub %in% to_remove,]
write.csv(share_cuticle_can, "candidates_sharing_eQTLs.csv", row.names=F)



### prepare the datasets needed for genomic links
trait_pos <- trait[c(match(share_cuticle_can$query, trait$gene), match(share_cuticle_can$sub, trait$gene)),]; rownames(trait_pos) <- 1:nrow(trait_pos)
trait_pos$gene_chr <- paste0("Chr ", trait_pos$gene_chr)
trait_pos_uniq <- trait_pos[!duplicated(trait_pos),]
write.csv(trait_pos_uniq, "candidate_sharing_eQTLs_traitInfo.csv", row.names=F)

trait_pos_uniq <- read.csv("candidate_sharing_eQTLs_traitInfo.csv")


eQTL_pos <- rbind(setNames(share_cuticle_can[, c("chr", "query_start", "query_end")], c("qtl_chr", "qtl_start", "qtl_end")),
				  setNames(share_cuticle_can[, c("chr", "sub_start", "sub_end")], c("qtl_chr", "qtl_start", "qtl_end")))
rownames(eQTL_pos) <- 1:nrow(eQTL_pos)
eQTL_pos$qtl_chr <- paste0("Chr ", eQTL_pos$qtl_chr)

#### use rgb() to generate the color with desired transparency and col2rgb() to find the code for desired colors, e.g., col2rgb("lightblue")
### color for cuticle candidates
col1 <- rgb(col2rgb("red")[1,], col2rgb("red")[2,], col2rgb("red")[3,], max = 255, alpha = 50)
### color for the other candidates
col2 <- rgb(col2rgb("blue")[1,], col2rgb("blue")[2,], col2rgb("blue")[3,], max = 255, alpha = 250)

pdf("trait_eQTL_shared_cuticle_candidates.pdf", width=6.8, height=6.8, pointsize=10)
circos.genomicInitialize(genome, axis.labels.cex=0.8, labels.cex=1)

circos.genomicTrack(trait_pos_uniq[!trait_pos_uniq$gene %in% cuticle_can, 4:7], stack=F, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicPoints(region, value, pch=16, col="red", cex=1,...)
})

labels <- trait_pos_uniq[!trait_pos_uniq$gene %in% cuticle_can, c(4:6, 3)]
#labels <- trait_pos[99:nrow(trait_pos), c(2:4, 1)]; 
#labels <- labels[!labels$gene %in%cuticle_can,]
#labels <- labels[!duplicated(labels),]; labels[,4] <- gsub("Zm00001d0", "", labels[,4])

circos.genomicLabels(labels, labels.column=4, side="inside", cex=0.8, connection_height = mm_h(2))

circos.genomicTrack(trait_pos_uniq[trait_pos_uniq$gene %in% cuticle_can, 4:7], stack=F, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicPoints(region, value, pch=16, col="blue", cex=1,...)
})
#labels <- trait_pos[1:98, c(2:4, 1)]; labels <- labels[!duplicated(labels),]; labels[,4] <- gsub("Zm00001d0", "", labels[,4])
labels <- trait_pos_uniq[trait_pos_uniq$gene %in% cuticle_can, c(4:6, 3)]
circos.genomicLabels(labels, labels.column=4, side="inside", cex=0.8, connection_height = mm_h(1))

circos.genomicTrack(snp, stack=T, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicText(region, value, labels="I", cex=1,...)
})
circos.genomicLink(trait_pos[,2:4], eQTL_pos, col=c(rep(col2, nrow(eQTL_pos)/2), rep(col1, nrow(eQTL_pos)/2)), lwd=0.2, border=c(rep("blue", nrow(eQTL_pos)/2), rep("red", nrow(eQTL_pos)/2)), direction=-1, arr.length=0.3)

dev.off()

#### create the table for qtls shared by cuticle candidates and the other candidates
share_cuticle_can1 <- share_cuticle_can[,c("query", "chr", "query_start", "query_end", "sub", "sub_start", "sub_end")]
colnames(share_cuticle_can1)[2] <- "qtl_chr"
query_annotation <- sub_annotation <- query_name <- sub_name <- NULL
query_start <- query_end <- sub_start <- sub_end <- NULL
query_chr <- sub_chr <- NULL
query_peak_lod <- sub_peak_lod <- NULL

for(i in 1:nrow(share_cuticle_can1)){
	query <- as.vector(share_cuticle_can1[i,"query"])
	sub <- as.vector(share_cuticle_can1[i,"sub"])
	
	query_annotation <- c(query_annotation, as.vector(trait_pos_uniq[trait_pos_uniq$gene == query, "annotation"]))
	sub_annotation <- c(sub_annotation, as.vector(trait_pos_uniq[trait_pos_uniq$gene == sub, "annotation"]))
	query_name <- c(query_name, as.vector(trait_pos_uniq[trait_pos_uniq$gene == query, "name"]))
	sub_name <- c(sub_name, as.vector(trait_pos_uniq[trait_pos_uniq$gene == sub, "name"]))
	
	query_start <- c(query_start, trait_pos_uniq[trait_pos_uniq$gene == query, "gene_start"])
	query_end <- c(query_end, trait_pos_uniq[trait_pos_uniq$gene == query, "gene_end"])
	sub_start <- c(sub_start, trait_pos_uniq[trait_pos_uniq$gene == sub, "gene_start"])
	sub_end <- c(sub_end, trait_pos_uniq[trait_pos_uniq$gene == sub, "gene_end"])
	query_chr <- c(query_chr, trait_pos_uniq[trait_pos_uniq$gene == query, "gene_chr"])
	sub_chr <- c(sub_chr, trait_pos_uniq[trait_pos_uniq$gene == sub, "gene_chr"])
	
	query_peak_lod <- c(query_peak_lod, trait_pos_uniq[trait_pos_uniq$gene == query, "peak_lod"])
	sub_peak_lod <- c(sub_peak_lod, trait_pos_uniq[trait_pos_uniq$gene == sub, "peak_lod"])		
}

share_cuticle_can1 <- data.frame(query=share_cuticle_can1$query, query_annotation, query_name, query_chr, query_start, query_end,
           sub=share_cuticle_can1$sub, sub_annotation, sub_name, sub_chr, sub_start, sub_end,
           qtl_chr=share_cuticle_can1[,"qtl_chr"],
           query_qtl_start=share_cuticle_can1[, "query_start"], query_qtl_end=share_cuticle_can1[, "query_end"], query_peak_lod,
           sub_qtl_start=share_cuticle_can1[, "sub_start"], sub_qtl_end=share_cuticle_can1[, "sub_end"], sub_peak_lod)
write.csv(share_cuticle_can1, "candidates_sharing_eQTLs.csv", row.names=F)


#### eQTLs for cuticle candidates only
pdf("trait_eQTL_cuticle_candidates.pdf", width=6.8, height=6.8, pointsize=10)
col2 <- rgb(col2rgb("blue")[1,], col2rgb("blue")[2,], col2rgb("blue")[3,], max = 255, alpha = 250)
circos.genomicInitialize(genome, axis.labels.cex=0.8, labels.cex=1)

circos.genomicTrack(trait_pos_uniq[trait_pos_uniq$gene %in% cuticle_can, 4:7], stack=F, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicPoints(region, value, pch=16, col="blue", cex=1,...)
})
labels <- trait_pos_uniq[trait_pos_uniq$gene %in% cuticle_can, c(4:6, 3)]
circos.genomicLabels(labels, labels.column=4, side="inside", cex=0.8, connection_height = mm_h(2))
circos.genomicTrack(snp, stack=T, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicText(region, value, labels="I", cex=1,...)
})
circos.genomicLink(trait_pos[trait_pos$gene %in% cuticle_can, 2:4], eQTL_pos[trait_pos$gene %in% cuticle_can, ], col=col2, lwd=0.2, border="blue", direction=-1, arr.length=0.3)
dev.off()






#### finding candidates with overlapping eQTLs with the cuticle non-candidates
share_cuticle_noncan <- NULL
cuticle_non_can <- cuticle_eQTL[!cuticle_eQTL$gene %in% unique(cuticle_can_eQTL[,"gene"]),]
for(i in 1:nrow(cuticle_non_can)){
	query <- cuticle_non_can[i, "gene"]
	query_chr <- cuticle_non_can[i, "qtl_chr"]
	query_start <- cuticle_non_can[i, "qtl_start"]
	query_end <- cuticle_non_can[i, "qtl_end"]
	share_cuticle_noncan <-  rbind(share_cuticle_noncan, share_eQTL(query=query, query_chr=query_chr, query_start=query_start, query_end=query_end, subject=candidate_eQTL))
}
to_remove <- c("Zm00001d018309", "Zm00001d048311")
share_cuticle_noncan <- share_cuticle_noncan[!share_cuticle_noncan$sub %in% to_remove, ]

### prepare the datasets needed for genomic links
trait_pos <- trait[c(match(share_cuticle_noncan$query, trait$gene), match(share_cuticle_noncan$sub, trait$gene)),]; rownames(trait_pos) <- 1:nrow(trait_pos)
trait_pos$gene_chr <- paste0("Chr", trait_pos$gene_chr)
eQTL_pos <- rbind(setNames(share_cuticle_noncan[, c("chr", "query_start", "query_end")], c("qtl_chr", "qtl_start", "qtl_end")),
				  setNames(share_cuticle_noncan[, c("chr", "sub_start", "sub_end")], c("qtl_chr", "qtl_start", "qtl_end")))
rownames(eQTL_pos) <- 1:nrow(eQTL_pos)
eQTL_pos$qtl_chr <- paste0("Chr", eQTL_pos$qtl_chr)


### color for cuticle lipids (not candidates)
col1 <- rgb(col2rgb("yellow")[1,], col2rgb("yellow")[2,], col2rgb("yellow")[3,], max = 255, alpha = 200)
### color for the other candidates
col2 <- rgb(col2rgb("blue")[1,], col2rgb("blue")[2,], col2rgb("blue")[3,], max = 255, alpha = 10)

pdf("trait_eQTL_lipid_orthologs.pdf", width=7.5, height=7.5)
circos.genomicInitialize(genome, axis.labels.cex=0.6, labels.cex=1)
circos.genomicTrack(snp, stack=T, track.height=0.05, panel.fun=function(region, value, ...){
				circos.genomicText(region, value, labels="I", cex=1,...)
})
circos.genomicTrack(trait_pos[1:443, 2:4], stack=T, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicPoints(region, value, pch=16, col="orange", cex=0.5,...)
})
circos.genomicTrack(trait_pos[444:nrow(trait_pos), 2:4], stack=T, track.height=0.05, panel.fun=function(region, value, ...){
		circos.genomicPoints(region, value, pch=16, col="blue", cex=0.5,...)
})
circos.genomicLink(trait_pos[,2:4], eQTL_pos, col=c(rep(col1, 443), rep(col2, 443)), lwd=0.2, border=c(rep("orange", 443), rep("blue", 443))
dev.off()
circos.genomicLink(trait_pos[1:443,2:4], eQTL_pos[1:443,], col=c(rep(col1, 443)), lwd=0.2, border=c(rep("orange", 443)))
