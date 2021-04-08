library(ggplot2)
library(gplots)
library(dplyr)
library(data.table)
library(vdjbaseVis)
library(fastmatch)
library(reshape2)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_3/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")
sources_folder <- paste0(project_folder, "figures/sources/")
  
# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

# load sources
source(paste0(sources_folder, "genotype_heatmap/internal_functions.R"))
source(paste0(sources_folder, "genotype_heatmap/genoHeatmap.R"))

# load DS1 genotypes after filtering unreliable unknown alleles
genotypes <- read.delim(paste0(required_files_folder, "HCV_Genotypes_filtered.tab"), sep = "\t", stringsAsFactors = F)

# Change the names of the collapsed genes
# TRBV6-23 > TRBV6-2/TRBV6-3
GENE.loc[["TRB"]] <- gsub("TRBV6-23", "TRBV6-2/6-3", GENE.loc[["TRB"]])
genotypes$GENE <- gsub("TRBV6-23", "TRBV6-2/6-3", genotypes$GENE)

##########################################################################################################################
############################ Genrate genotype heatmaps ###################################################################
##########################################################################################################################

genotypes$SUBJECT[grepl("^C[0-9]", genotypes$SUBJECT)] <- paste0("H", genotypes$SUBJECT[grepl("^C[0-9]", genotypes$SUBJECT)])
genotypes$SUBJECT <- gsub("T", "", genotypes$SUBJECT)

genotypes <- genotypes[!genotypes$GENE %in% PSEUDO[["TRB"]],]

genotypes$GENOTYPED_ALLELES <- gsub("Deletion", "Del", genotypes$GENOTYPED_ALLELES)
alleles <- unique(unlist(strsplit(genotypes$GENOTYPED_ALLELES, ",", fixed = T)))
genotypes$GENOTYPED_ALLELES <- gsub("_[0-9]+[A-Z]+[0-9]+", "", genotypes$GENOTYPED_ALLELES)

genotypes$ALLELE_FACTOR <- lapply(1:nrow(genotypes), function(x) {rep(0, length(alleles))})
for (i in 1:length(alleles)) {
  allele <- alleles[[i]]
  allele_factor <- rep(0, length(alleles))
  allele_factor[i] <- 1
  genotypes$ALLELE_FACTOR[grepl(paste0(allele,"(,|$)"), genotypes$GENOTYPED_ALLELES)] <- lapply(
    genotypes$ALLELE_FACTOR[grepl(paste0(allele,"(,|$)"), genotypes$GENOTYPED_ALLELES)], function(x){
      x + allele_factor
    })
}

genotypes_matrix <- genotypes %>% select(SUBJECT, GENE, ALLELE_FACTOR)
genotypes_matrix <- dcast(genotypes_matrix, SUBJECT ~ GENE)

genotypes_matrix[is.na(genotypes_matrix)] <- lapply(genotypes_matrix[is.na(genotypes_matrix)], function(x) {rep(0, length(alleles))})
subject_distances_df <- data_frame(SUBJECT=character(), COMPARE_TO=character(), DISTANCE=numeric())
for (subject1 in genotypes_matrix$SUBJECT) {
  for (subject2 in genotypes_matrix$SUBJECT) {
    distance <- 0
    for (col in 2:ncol(genotypes_matrix)) {
      gene_dis <- abs(unlist(genotypes_matrix[genotypes_matrix$SUBJECT == subject1, col]) - 
                        unlist(genotypes_matrix[genotypes_matrix$SUBJECT == subject2, col]))
      gene_dis <- sum(gene_dis)
      distance <- distance + gene_dis
    }
    subject_distances_df[nrow(subject_distances_df) + 1,] <- list(subject1, subject2, distance)
  }
}


subject_distances_df$DISTANCE <- as.numeric(subject_distances_df$DISTANCE)
subject_distances <- dcast(subject_distances_df, SUBJECT ~ COMPARE_TO)
rownames(subject_distances) <- subject_distances$SUBJECT 
subject_distances[,-1] <- as.numeric(subject_distances[,-1])
subject_distances <- as.matrix(subject_distances[,-1])

hcr <- hclust(dist(subject_distances))
ddr <- as.dendrogram(hcr)
hcr$order

Rowv <- rowMeans(subject_distances, na.rm = T)
ddr <- reorder(ddr, Rowv)

subject_order <- rownames(subject_distances)[order.dendrogram(ddr)]
genotypes$SUBJECT <- factor(genotypes$SUBJECT, levels=subject_order)
genotypes <- genotypes[order(genotypes$SUBJECT),]


v_genotypes <- genotypes[grepl("TRBV",genotypes$GENE),]
v_genotype_heatmap <- genoHeatmap(geno_table = v_genotypes, chain = "TRB", lk_cutoff = 3, pseudo_genes = T, ORF_genes = F, removeIGH = F)

pdf(file = paste0(figure_folder, "TRBV_DS1_genotype.pdf"), height = 10, width = v_genotype_heatmap$width)
print(v_genotype_heatmap$p)
dev.off()

# png(file = paste0(figure_folder, "TRBV_DS1_genotype.png"))
# print(v_genotype_heatmap$p)
# dev.off()

dj_genotypes <- genotypes[grepl("TRBD|J",genotypes$GENE),]
dj_genotype_heatmap <- genoHeatmap(geno_table = dj_genotypes, chain = "TRB", lk_cutoff = 3, pseudo_genes = F, removeIGH = F)


pdf(file = paste0(figure_folder, "TRBDJ_DS1_genotype.pdf"), height = 10, width = dj_genotype_heatmap$width)
print(dj_genotype_heatmap$p)
dev.off()

# png(file = paste0(figure_folder, "TRBDJ_DS1_genotype.png"))
# print(dj_genotype_heatmap$p)
# dev.off()

