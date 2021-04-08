library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(cowplot)
library(magrittr)
library(gridExtra)
library(ggpubr)

##################################################################################################################################
############################ Load required files #################################################################################
##################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_14/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

# read genotypes 
genotypes <- read.delim(paste0(required_files_folder, "BIOMED2_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)
v_frequencies <- read.delim(paste0(required_files_folder, "BIOMED2_gene_usage.tab"), header = T, sep = "\t", stringsAsFactors = F)

##################################################################################################################################################
###########################  Figure  #############################################################################################################
##################################################################################################################################################

# check allele distribution
homozygous <- genotypes[!grepl(",", genotypes$GENOTYPED_ALLELES),]

homozygous <- homozygous %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N = n()*2)

heterozygous <- genotypes[grepl(",", genotypes$GENOTYPED_ALLELES),]
heterozygous <- heterozygous %>% separate_rows(GENOTYPED_ALLELES, sep = ",")
heterozygous <- heterozygous %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N = n())

allele_dis <- rbind(homozygous, heterozygous)
allele_dis <- allele_dis %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N = sum(N))
allele_dis$TOTAL_GENE <- unlist(lapply(allele_dis$GENE, function(gene) {sum(allele_dis$N[allele_dis$GENE == gene])}))
allele_dis$FREQ <- allele_dis$N / allele_dis$TOTAL_GENE

# check allelic realations between genes
all_genes <- unique(genotypes$GENE)
genotypes$GENOTYPED_ALLELES[grepl("Del", genotypes$GENOTYPED_ALLELES)] <- "Del"
genes2check <- unlist(lapply(all_genes,function(gene){if(max(allele_dis$FREQ[allele_dis$GENE==gene]) < 0.7){gene}}))
genes2check <- genes2check[!grepl("OR",genes2check)]
# present the genotype in spesific format

genotypes$GENOTYPE <- unlist(lapply(genotypes$GENOTYPED_ALLELES, function(x) {
  x <- unlist(strsplit(x, ",", fixed = T));
  x <- x[order(x)];
  paste(x, collapse = "-")
}))
gene<-"TRBV7-2"


v_frequencies <- v_frequencies[!grepl("[a-z]$", v_frequencies$SUBJECT),]

v_frequencies$GENE <- gsub("TRB", "", v_frequencies$GENE)
TRBV_LOC <- GENE.loc[["TRB"]][grepl("V", GENE.loc[["TRB"]])]
TRBV_LOC <- TRBV_LOC[!grepl("5-2|7-5|22-1|V26|V8", TRBV_LOC)]

TRBV_LOC <- TRBV_LOC[grepl("V5-1|V6-1|V4-2|6-23|V3-2|V4-3|V7-2", TRBV_LOC)]

GENE_GENO <- genotypes[genotypes$GENE == "TRBV7-2" & genotypes$K_DIFF >= 3,]
GENE_GENO <- GENE_GENO[!grepl("[a-d]$", GENE_GENO$SUBJECT),]

GENE_GENO <- GENE_GENO %>% select(SUBJECT, GENOTYPE)

v_frequencies <- merge(v_frequencies, GENE_GENO, by="SUBJECT")
v_frequencies <- v_frequencies[!grepl("A251G|C135G", v_frequencies$GENOTYPE),]
v_frequencies$GENOTYPE <- as.factor(v_frequencies$GENOTYPE)

v_frequencies <- v_frequencies[v_frequencies$GENE %in% gsub("TRB", "", TRBV_LOC),]
v_frequencies$GENE <- factor(v_frequencies$GENE, levels = gsub("TRB", "", TRBV_LOC))


v_usage_graph <- ggplot(v_frequencies, aes(x=GENE, y=FREQ, fill = GENOTYPE)) + 
  geom_boxplot() +  
  # scale_x_discrete(limits=gsub("TRB", "", TRBV_LOC)) + 
  ylab("Usage") + xlab("TRBV gene") +
  guides(fill=guide_legend("TRBV7-2 genotype", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14), axis.title = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        legend.position=c(0.5,0.8))


comperasiments <- list(c("bp01", "bp01-bp02"), c("bp02", "bp01-bp02"), c("bp01", "bp02"))
xmin_space <- c(-0.3, 0.025, -0.3)
xmax_space <- c(-0.025, 0.3, 0.3)
y.values <- c(0.205, 0.205, 0.22)

i <- 1
for (i in 1:3) {
  temp <- v_frequencies[v_frequencies$GENOTYPE %in% comperasiments[[i]],] 
  p.values <- sapply(split(temp, temp$GENE), function(x){wilcox.test(FREQ~GENOTYPE, x)$p.value})
  p.values <- p.adjust(p.values, method = "bonferroni")
  
  labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*",paste0(ifelse(i==2, " ", ""), "n.s.", ifelse(i==1, " ", ""))))
  # if (i == 3) {
  #   y.values <- y.values + 0.025
  # }
  v_usage_graph <- v_usage_graph + geom_signif(y_position = y.values[i], xmin = 1:length(TRBV_LOC) + xmin_space[i], xmax = 1:length(TRBV_LOC) + xmax_space[i], annotations = labels, textsize = 4)
}

ggsave(paste0(figure_folder, "TRBV7-2_partitial.pdf"), v_usage_graph, height = 8, width = 14)
# ggsave(paste0(figure_folder, "TRBV7-2_partitial.png"), v_usage_graph, height = 8, width = 14)


