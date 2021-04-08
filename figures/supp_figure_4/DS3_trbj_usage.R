library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_4/")

DS3_genotypes <- read.delim(paste0(required_files_folder, "BIOMED2_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)

#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

DS3_genotypes <- DS3_genotypes[!grepl("[a-d]$", DS3_genotypes$SUBJECT),]

ALL_J_GENO <- DS3_genotypes[grepl("TRBJ", DS3_genotypes$GENE),]

ALL_GENO_TRBJ1_6 <- DS3_genotypes[DS3_genotypes$GENE=="TRBJ1-6",]
trbj1_6_01_subjects <- ALL_GENO_TRBJ1_6$SUBJECT[ALL_GENO_TRBJ1_6$GENOTYPED_ALLELES == "01"]
trbj1_6_02_subjects <- ALL_GENO_TRBJ1_6$SUBJECT[ALL_GENO_TRBJ1_6$GENOTYPED_ALLELES == "02"]
trbj1_6_hetero_subjects <- ALL_GENO_TRBJ1_6$SUBJECT[grepl(",", ALL_GENO_TRBJ1_6$GENOTYPED_ALLELES, fixed = T)]


ALL_J_GENO$J1_6_GENO <- "J1_6_01_02"
ALL_J_GENO$J1_6_GENO[ALL_J_GENO$SUBJECT %in% trbj1_6_01_subjects] <- "J1_6_01"
ALL_J_GENO$J1_6_GENO[ALL_J_GENO$SUBJECT %in% trbj1_6_02_subjects] <- "J1_6_02"


total_trbj <- aggregate(ALL_J_GENO$TOTAL, by=list(SUBJECT=ALL_J_GENO$SUBJECT), FUN=sum)
names(total_trbj)[names(total_trbj) == "x"] <- "TOTAL_J_SEQS"
ALL_J_GENO <- merge(ALL_J_GENO, total_trbj, by = "SUBJECT")

ALL_J_GENO$J_FREQ <- ALL_J_GENO$TOTAL / ALL_J_GENO$TOTAL_J_SEQS

ALL_J_GENO$GENE <- gsub("TRB", "", ALL_J_GENO$GENE)
biomed_trbj_usage <- ggplot(ALL_J_GENO, aes(x=GENE, y=J_FREQ, fill = J1_6_GENO)) + 
  geom_boxplot() + 
  ylab("Usage") + xlab("TRBJ gene") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14), axis.title = element_text(size=18))


# ggsave(paste0(figure_folder, "BIOMED_TRBJ_USAGE.png"), biomed_trbj_usage, height = 6)
ggsave(paste0(figure_folder, "BIOMED_TRBJ_USAGE.pdf"), biomed_trbj_usage, width = 20, height = 6)
