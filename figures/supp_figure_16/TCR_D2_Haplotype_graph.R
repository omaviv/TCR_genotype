library(seqinr)
library(tigger)
library(optparse)
library(gtools)
library(ggplot2)
library(utils)
library(stringr)
library(dplyr)
library(data.table)
library(alakazam)
library(reshape2)
library(rabhit)
library(cowplot)
library(fastmatch)

##################################################################################################################################################
############################ Load required files #################################################################################################
##################################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_16/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")
sources_folder <- paste0(project_folder, "figures/sources/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

# load sources
source(paste0(sources_folder, "haplotype/createFullHaplotype.R"))
source(paste0(sources_folder, "haplotype/internal_functions.R"))
source(paste0(sources_folder, "haplotype/Haplotype_graph.R"))


DS1_haplotypes <- read.delim(paste0(required_files_folder, "HCV_Haplotypes_D2_full.tab"), sep = "\t", stringsAsFactors = F)

# Change the names of the collapsed genes
# TRBV6-23 > TRBV6-2/TRBV6-3
GENE.loc[["TRB"]] <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", GENE.loc[["TRB"]])
DS1_haplotypes$gene <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", DS1_haplotypes$gene)

###################################################################################################################################
############################ Generate haplotype graph #############################################################################
###################################################################################################################################

names(DS1_haplotypes) <- toupper(names(DS1_haplotypes))
DS1_haplotypes$SUBJECT[grepl("^C[0-9]", DS1_haplotypes$SUBJECT)] <- paste0("H", DS1_haplotypes$SUBJECT[grepl("^C[0-9]", DS1_haplotypes$SUBJECT)])
DS1_haplotypes$SUBJECT <- gsub("T", "", DS1_haplotypes$SUBJECT)

haplo_d2 <- DS1_haplotypes
haplo_d2 <- haplo_d2[!haplo_d2$GENE %in% PSEUDO[["TRB"]],]
haplo_d2 <- haplo_d2[!grepl("OR", haplo_d2$GENE),]

haplo_d2$TRBD2_01 <-  gsub("_[0-9]+[A-Z]+[0-9]+", "", haplo_d2$TRBD2_01)
haplo_d2$TRBD2_02 <-  gsub("_[0-9]+[A-Z]+[0-9]+", "", haplo_d2$TRBD2_02)
haplo_d2$ALLELES <-  gsub("_[0-9]+[A-Z]+[0-9]+", "", haplo_d2$ALLELES)

haplo_d2_graph <- hapHeatmap(hap_table = haplo_d2, chain = "TRB", removeIGH = F, lk_cutoff = 3)

pdf(file = paste0(figure_folder, "DS1_Haplotype_D2.pdf"), height = 18, width = haplo_d2_graph$width)
print(haplo_d2_graph$p)
dev.off()

# png(file = paste0(figure_folder, "DS1_DS2_Haplotype_D2.pdf"))
# print(haplo_d2_graph$p)
# dev.off()

