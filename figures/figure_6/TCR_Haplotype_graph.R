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
figure_folder <- paste0(project_folder, "figures/figure_6/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")
sources_folder <- paste0(project_folder, "figures/sources/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

# load sources
source(paste0(sources_folder, "haplotype/createFullHaplotype.R"))
source(paste0(sources_folder, "haplotype/internal_functions.R"))
source(paste0(sources_folder, "haplotype/Haplotype_graph.R"))


DS1_haplotypes <- read.delim(paste0(required_files_folder, "HCV_Haplotypes_J1_6_full.tab"), sep = "\t", stringsAsFactors = F)
DS2_haplotypes <- read.delim(paste0(required_files_folder, "SC_Haplotypes_J1_6_full.tab"), sep = "\t", stringsAsFactors = F)

# Change the names of the collapsed genes
# TRBV6-23 > TRBV6-2/TRBV6-3
GENE.loc[["TRB"]] <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", GENE.loc[["TRB"]])
DS1_haplotypes$gene <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", DS1_haplotypes$gene)
DS2_haplotypes$gene <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", DS2_haplotypes$gene)

###################################################################################################################################
############################ Generate haplotype graph #############################################################################
###################################################################################################################################

names(DS1_haplotypes) <- toupper(names(DS1_haplotypes))
DS1_haplotypes$SUBJECT[grepl("^C[0-9]", DS1_haplotypes$SUBJECT)] <- paste0("H", DS1_haplotypes$SUBJECT[grepl("^C[0-9]", DS1_haplotypes$SUBJECT)])
DS1_haplotypes$SUBJECT <- gsub("T", "", DS1_haplotypes$SUBJECT)
DS1_haplotypes$SUBJECT <- paste0("DS1:", DS1_haplotypes$SUBJECT)

names(DS2_haplotypes) <- toupper(names(DS2_haplotypes))
DS2_haplotypes <- DS2_haplotypes[grepl("TRB", DS2_haplotypes$GENE),]
DS2_haplotypes$SUBJECT <- paste0("DS2:", DS2_haplotypes$SUBJECT)

haplo_j1_6 <- rbind(DS1_haplotypes, DS2_haplotypes)
haplo_j1_6 <- haplo_j1_6[!haplo_j1_6$GENE %in% PSEUDO[["TRB"]],]

names(haplo_j1_6) <- gsub(".", "-", names(haplo_j1_6), fixed = T)
haplo_j1_6$`TRBJ1-6_01` <-  gsub("_[0-9]+[A-Z]+[0-9]+", "", haplo_j1_6$`TRBJ1-6_01`)
haplo_j1_6$`TRBJ1-6_02` <-  gsub("_[0-9]+[A-Z]+[0-9]+", "", haplo_j1_6$`TRBJ1-6_02`)

haplo_j1_6_graph <- hapHeatmap(hap_table = haplo_j1_6, chain = "TRB", removeIGH = F, lk_cutoff = 3)

pdf(file = paste0(figure_folder, "DS1_DS2_Haplotype_J1-6.pdf"), height = 15, width = haplo_j1_6_graph$width)
print(haplo_j1_6_graph$p)
dev.off()

# png(file = paste0(figure_folder, "DS1_DS2_Haplotype_J1-6.pdf"))
# print(haplo_j1_6_graph$p)
# dev.off()

