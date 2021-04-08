library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(wesanderson)


##################################################################################################################################
############################ Load required files #################################################################################
##################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_13/")

genotypes <- read.delim(paste0(required_files_folder, "BIOMED2_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)

##################################################################################################################################################
###########################  Figure  #############################################################################################################
##################################################################################################################################################

# read the merged genotype file
genotypes <- genotypes[!grepl("[a-d]$", genotypes$SUBJECT),]
genotypes <- genotypes[genotypes$K_DIFF >= 3,]
genotypes$GENOTYPED_ALLELES[grepl("Del", genotypes$GENOTYPED_ALLELES)] <- "Del"

# present the genotype in spesific format
genotypes$GENOTYPE <- unlist(lapply(genotypes$GENOTYPED_ALLELES, function(x) {
  x <- unlist(strsplit(x, ",", fixed = T));
  x <- x[order(x)];
  paste(x, collapse = "-")
}))

allele_relations <- genotypes %>% select(GENE, GENOTYPE, SUBJECT)
allele_relations <- merge(allele_relations, allele_relations, by = "SUBJECT")
allele_relations <- allele_relations %>% group_by(GENE.x, GENE.y, GENOTYPE.x, GENOTYPE.y) %>% dplyr::summarise(N=n())
genes <- unique(allele_relations$GENE.x)

# calculate the frequency
allele_relations$FREQ <- 0
for (i in 1:(length(genes))) {
  gene1 <- genes[i]
  for (j in i:length(genes)) {
    gene2 <- genes[j]
    total <- sum(allele_relations$N[allele_relations$GENE.x == gene1 & allele_relations$GENE.y == gene2])
    allele_relations$FREQ[allele_relations$GENE.x == gene1 & allele_relations$GENE.y == gene2] <- allele_relations$N[allele_relations$GENE.x == gene1 & allele_relations$GENE.y == gene2] / total
    allele_relations$FREQ[allele_relations$GENE.x == gene2 & allele_relations$GENE.y == gene1] <- allele_relations$N[allele_relations$GENE.x == gene2 & allele_relations$GENE.y == gene1] / total
  }
}

# remove genes that the maximum square frequency of one of its alleles is above the threshold 
max_allele_homozygous_dis_threshold <- 0.9
genes2remove <- c()
same_genes <- allele_relations[allele_relations$GENE.x == allele_relations$GENE.y,]
for (gene in genes) {
  temp <- same_genes[same_genes$GENE.x == gene,]
  if (max(temp$FREQ) >= max_allele_homozygous_dis_threshold) {
    genes2remove <- c(genes2remove, gene)
  }
}

allele_relations <- allele_relations[!(allele_relations$GENE.x %in% genes2remove) &
                                       !(allele_relations$GENE.y %in% genes2remove),]
allele_relations <- allele_relations[allele_relations$GENE.x != allele_relations$GENE.y,]

allele_relations$LABELS <- paste0(allele_relations$N, "\n(", round(allele_relations$FREQ, 4), ")")

gene2 <- "TRBV7-2"
gene1 <- "TRBV4-3"
temp <- allele_relations[allele_relations$GENE.x==gene1 & allele_relations$GENE.y==gene2,]

allele_relation_graph <- ggplot(temp, aes(x = GENOTYPE.x, y = GENOTYPE.y, fill = FREQ)) +
  geom_tile() + scale_fill_continuous(limits=c(0, 1), low = "white", high = "darkgrey") +
  geom_text(aes(label = LABELS)) +
  xlab("TRBV4-3 genotype") + ylab("TRBV7-2 genotype") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=8), plot.title = element_text(hjust = 0.5))


ggsave(paste0(figure_folder, "BIOMED_allele_relations.pdf"), allele_relation_graph, width = 10, height = 9, limitsize = F)
