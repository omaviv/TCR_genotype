library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
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
figure_folder <- paste0(project_folder, "figures/supp_figure_9/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

DS3_genotypes <- read.delim(paste0(required_files_folder, "BIOMED2_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)
DS4_genotypes <- read.delim(paste0(required_files_folder, "Adaptive_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)

##################################################################################################################################
############################ Calculate D2 genotype boundaries ####################################################################
##################################################################################################################################
# general boundaries
Adaptive_GENO_TRBD2 <- DS4_genotypes[DS4_genotypes$GENE=="TRBD2",]

Adaptive_GENO_TRBD2 <- separate(Adaptive_GENO_TRBD2, "ALLELES", paste("ALLELES", 1:2, sep = "_"), sep = ",")
Adaptive_GENO_TRBD2 <- separate(Adaptive_GENO_TRBD2, "COUNTS", paste("COUNTS", 1:2, sep = "_"), sep = ",")

Adaptive_GENO_TRBD2$COUNT_01 <- 0
Adaptive_GENO_TRBD2$COUNT_01[grepl("01", Adaptive_GENO_TRBD2$ALLELES_1)] <- Adaptive_GENO_TRBD2$COUNTS_1[grepl("01", Adaptive_GENO_TRBD2$ALLELES_1)] 
Adaptive_GENO_TRBD2$COUNT_01[grepl("01", Adaptive_GENO_TRBD2$ALLELES_2)] <- Adaptive_GENO_TRBD2$COUNTS_2[grepl("01", Adaptive_GENO_TRBD2$ALLELES_2)] 

Adaptive_GENO_TRBD2$COUNT_01 <- as.numeric(Adaptive_GENO_TRBD2$COUNT_01)
Adaptive_GENO_TRBD2$FREQ_01 <- Adaptive_GENO_TRBD2$COUNT_01 / Adaptive_GENO_TRBD2$TOTAL


D2_01_HOMO_FREQ <- Adaptive_GENO_TRBD2$FREQ_01[Adaptive_GENO_TRBD2$FREQ_01 > 0.75]
D2_HETERO_FREQ <- Adaptive_GENO_TRBD2$FREQ_01[Adaptive_GENO_TRBD2$FREQ_01 > 0.25 & Adaptive_GENO_TRBD2$FREQ_01 < 0.75]
D2_02_HOMO_FREQ <- Adaptive_GENO_TRBD2$FREQ_01[Adaptive_GENO_TRBD2$FREQ_01 < 0.25]

D2_01_HOMO_HETERO_BOUND <- (mean(D2_01_HOMO_FREQ)*sd(D2_HETERO_FREQ) + mean(D2_HETERO_FREQ)*sd(D2_01_HOMO_FREQ)) / (sd(D2_01_HOMO_FREQ)+sd(D2_HETERO_FREQ))
D2_02_HOMO_HETERO_BOUND <- (mean(D2_02_HOMO_FREQ)*sd(D2_HETERO_FREQ) + mean(D2_HETERO_FREQ)*sd(D2_02_HOMO_FREQ)) / (sd(D2_02_HOMO_FREQ)+sd(D2_HETERO_FREQ))

# DS3 restricted boundaries
BIOMED2_GENO_TRBD2 <- DS3_genotypes[DS3_genotypes$GENE=="TRBD2",]
BIOMED2_GENO_TRBD2 <- BIOMED2_GENO_TRBD2[!grepl("[a-d]$", BIOMED2_GENO_TRBD2$SUBJECT),]

BIOMED2_GENO_TRBD2 <- separate(BIOMED2_GENO_TRBD2, "ALLELES", paste("ALLELES", 1:2, sep = "_"), sep = ",")
BIOMED2_GENO_TRBD2 <- separate(BIOMED2_GENO_TRBD2, "COUNTS", paste("COUNTS", 1:2, sep = "_"), sep = ",")

BIOMED2_GENO_TRBD2$COUNT_01 <- 0
BIOMED2_GENO_TRBD2$COUNT_01[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_1)] <- BIOMED2_GENO_TRBD2$COUNTS_1[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_1)] 
BIOMED2_GENO_TRBD2$COUNT_01[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_2)] <- BIOMED2_GENO_TRBD2$COUNTS_2[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_2)] 

BIOMED2_GENO_TRBD2$COUNT_01 <- as.numeric(BIOMED2_GENO_TRBD2$COUNT_01)
BIOMED2_GENO_TRBD2$FREQ_01 <- BIOMED2_GENO_TRBD2$COUNT_01 / BIOMED2_GENO_TRBD2$TOTAL

D2_HOMO02_BOUND <- (D2_02_HOMO_HETERO_BOUND + mean(D2_02_HOMO_FREQ)) / 2
D2_HETERO_1_BOUND <- (D2_02_HOMO_HETERO_BOUND + mean(D2_HETERO_FREQ)) / 2
D2_HETERO_2_BOUND <- (D2_01_HOMO_HETERO_BOUND + mean(D2_HETERO_FREQ)) / 2
D2_HOMO01_BOUND <- (D2_01_HOMO_HETERO_BOUND + mean(D2_01_HOMO_FREQ)) / 2

#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

BIOMED2_ALL_GENO <- DS3_genotypes[!grepl("[a-d]$", DS3_genotypes$SUBJECT),]
# BIOMED2_ALL_GENO <- BIOMED2_ALL_GENO[!grepl("hem", BIOMED2_ALL_GENO$SUBJECT),]
BIOMED2_GENO_TRBD2 <- BIOMED2_GENO_TRBD2[BIOMED2_GENO_TRBD2$FREQ_01 < D2_HOMO02_BOUND |
                                           (BIOMED2_GENO_TRBD2$FREQ_01 > D2_HETERO_1_BOUND & BIOMED2_GENO_TRBD2$FREQ_01 < D2_HETERO_2_BOUND)|
                                           BIOMED2_GENO_TRBD2$FREQ_01 > D2_HOMO01_BOUND,]

BIOMED2_GENO_TRBD2$D2_GENO <- "D2_01_02"
BIOMED2_GENO_TRBD2$D2_GENO[BIOMED2_GENO_TRBD2$FREQ_01 > D2_HOMO01_BOUND] <- "D2_01"
BIOMED2_GENO_TRBD2$D2_GENO[BIOMED2_GENO_TRBD2$FREQ_01 < D2_HOMO02_BOUND] <- "D2_02"


BIOMED2_GENO_TRBJ1_6 <- BIOMED2_ALL_GENO[BIOMED2_ALL_GENO$GENE == "TRBJ1-6",]
BIOMED2_GENO_TRBJ1_6 <- BIOMED2_GENO_TRBJ1_6[BIOMED2_GENO_TRBJ1_6$SUBJECT %in% BIOMED2_GENO_TRBD2$SUBJECT,]

BIOMED2_GENO_TRBJ1_6$J1_6_GENO <- "J1_6_01_02"
BIOMED2_GENO_TRBJ1_6$J1_6_GENO[BIOMED2_GENO_TRBJ1_6$GENOTYPED_ALLELES == "01"] <- "J1_6_01"
BIOMED2_GENO_TRBJ1_6$J1_6_GENO[BIOMED2_GENO_TRBJ1_6$GENOTYPED_ALLELES == "02"] <- "J1_6_02"

BIOMED2_GENO_TRBD2 <- BIOMED2_GENO_TRBD2 %>% select(SUBJECT, D2_GENO)
BIOMED2_GENO_TRBJ1_6 <- BIOMED2_GENO_TRBJ1_6 %>% select(SUBJECT, J1_6_GENO)

BIOMED2_TRBJ1_6_TRBD2 <- merge(BIOMED2_GENO_TRBD2, BIOMED2_GENO_TRBJ1_6, by = "SUBJECT")

BIOMED2_TRBJ1_6_TRBD2 <- BIOMED2_TRBJ1_6_TRBD2 %>% group_by(D2_GENO, J1_6_GENO) %>% dplyr::summarise(COUNT = n())

trbd2_geno <- c("D2_01", "D2_01_02", "D2_02")
trbj1_6_geno <- c("J1_6_01", "J1_6_01_02", "J1_6_02")

for (d2 in trbd2_geno) {
  for (j16 in trbj1_6_geno) {
    if (nrow(BIOMED2_TRBJ1_6_TRBD2[BIOMED2_TRBJ1_6_TRBD2$D2_GENO == d2 & BIOMED2_TRBJ1_6_TRBD2$J1_6_GENO == j16,]) == 0){
      BIOMED2_TRBJ1_6_TRBD2[nrow(BIOMED2_TRBJ1_6_TRBD2)+1,] <- list(d2, j16, 0)
    }
  }
}

BIOMED2_TRBJ1_6_TRBD2$FREQ <- BIOMED2_TRBJ1_6_TRBD2$COUNT / nrow(BIOMED2_GENO_TRBD2)

BIOMED2_TRBJ1_6_TRBD2$LABELS <- paste0(BIOMED2_TRBJ1_6_TRBD2$COUNT, "\n(", round(BIOMED2_TRBJ1_6_TRBD2$FREQ, 4), ")")

biomed2_j16_d2_relation_graph <- ggplot(BIOMED2_TRBJ1_6_TRBD2, aes(x = D2_GENO, y = J1_6_GENO, fill = FREQ)) +
  geom_tile() + scale_fill_continuous(limits=c(0, 1), low = "white", high = "darkgrey") +
  geom_text(aes(label = LABELS)) +
  xlab("TRBD2 Genotype") + ylab("TRBJ1-6 Genotype") + labs(fill="Intersection\nfrequency") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=8))

ggsave(paste0(figure_folder, "BIOMED2_D2_J1_6_RELATION.pdf"), biomed2_j16_d2_relation_graph)
# ggsave(paste0(figure_folder, "BIOMED2_D2_J1_6_RELATION.png"), biomed2_j16_d2_relation_graph)

