library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(magrittr)
library(gridExtra)
library(ggpubr)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_4/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

DS1_genotypes <- read.delim(paste0(required_files_folder, "HCV_Genotypes_const_2.tab"), header = T, sep = "\t", stringsAsFactors = F)
DS3_genotypes <- read.delim(paste0(required_files_folder, "BIOMED2_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)
DS4_genotypes <- read.delim(paste0(required_files_folder, "DS4_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)


##################################################################################################################################
############################ DS1 (HCV) TRBD2 ERROR GRAPH #########################################################################
##################################################################################################################################

HCV_GENO_TRBD2 <- DS1_genotypes[DS1_genotypes$GENE=="TRBD2",]

HCV_GENO_TRBD2 <- separate(HCV_GENO_TRBD2, "ALLELES", paste("ALLELES", 1:2, sep = "_"), sep = ",")
HCV_GENO_TRBD2 <- separate(HCV_GENO_TRBD2, "COUNTS", paste("COUNTS", 1:2, sep = "_"), sep = ",")

HCV_GENO_TRBD2$COUNT_01 <- 0
HCV_GENO_TRBD2$COUNT_01[grepl("01", HCV_GENO_TRBD2$ALLELES_1)] <- HCV_GENO_TRBD2$COUNTS_1[grepl("01", HCV_GENO_TRBD2$ALLELES_1)] 
HCV_GENO_TRBD2$COUNT_01[grepl("01", HCV_GENO_TRBD2$ALLELES_2)] <- HCV_GENO_TRBD2$COUNTS_2[grepl("01", HCV_GENO_TRBD2$ALLELES_2)] 

HCV_GENO_TRBD2$COUNT_01 <- as.numeric(HCV_GENO_TRBD2$COUNT_01)
HCV_GENO_TRBD2$FREQ_01 <- HCV_GENO_TRBD2$COUNT_01 / HCV_GENO_TRBD2$TOTAL

summary(HCV_GENO_TRBD2$TOTAL)

# d <- ggplot(HCV_GENO_TRBD2, aes(x=FREQ_01)) + 
#   geom_density(alpha=.2, fill="#FF6666") + ggtitle("TRBD2 Biomed2 dataset")

hcv_d2_error_graph <- ggplot(HCV_GENO_TRBD2, aes(x=FREQ_01)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.025), colour="black", fill="darkgrey") +
  xlab("TRBD2*01 fraction") + ylab("Number of individuals")+
  theme_classic() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=12), axis.text.x = element_text(size=12),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14))

hcv_d2_error_graph

require(grid)
title.grob <- textGrob(
  label = "A.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

hcv_d2_error_graph <- arrangeGrob(hcv_d2_error_graph, top = title.grob)

# ggsave(paste0(figure_folder, "DS1_TRBD2_allele_distribution.png"), hcv_d2_error_graph, height = 5)
# ggsave(paste0(figure_folder, "DS1_TRBD2_allele_distribution.pdf"), hcv_d2_error_graph, height = 5)

##################################################################################################################################
############################ DS4 (CMV) TRBD2 ERROR GRAPH #########################################################################
##################################################################################################################################

Adaptive_GENO_TRBD2 <- DS4_genotypes[DS4_genotypes$GENE=="TRBD2",]

Adaptive_GENO_TRBD2 <- separate(Adaptive_GENO_TRBD2, "ALLELES", paste("ALLELES", 1:2, sep = "_"), sep = ",")
Adaptive_GENO_TRBD2 <- separate(Adaptive_GENO_TRBD2, "COUNTS", paste("COUNTS", 1:2, sep = "_"), sep = ",")

Adaptive_GENO_TRBD2$COUNT_01 <- 0
Adaptive_GENO_TRBD2$COUNT_01[grepl("01", Adaptive_GENO_TRBD2$ALLELES_1)] <- Adaptive_GENO_TRBD2$COUNTS_1[grepl("01", Adaptive_GENO_TRBD2$ALLELES_1)] 
Adaptive_GENO_TRBD2$COUNT_01[grepl("01", Adaptive_GENO_TRBD2$ALLELES_2)] <- Adaptive_GENO_TRBD2$COUNTS_2[grepl("01", Adaptive_GENO_TRBD2$ALLELES_2)] 

Adaptive_GENO_TRBD2$COUNT_01 <- as.numeric(Adaptive_GENO_TRBD2$COUNT_01)
Adaptive_GENO_TRBD2$FREQ_01 <- Adaptive_GENO_TRBD2$COUNT_01 / Adaptive_GENO_TRBD2$TOTAL

summary(Adaptive_GENO_TRBD2$TOTAL)

# d <- ggplot(Adaptive_GENO_TRBD2, aes(x=FREQ_01)) + 
#   geom_density(alpha=.2, fill="#FF6666") + ggtitle("TRBD2 Biomed2 dataset")


D2_01_HOMO_FREQ <- Adaptive_GENO_TRBD2$FREQ_01[Adaptive_GENO_TRBD2$FREQ_01 > 0.75]
D2_HETERO_FREQ <- Adaptive_GENO_TRBD2$FREQ_01[Adaptive_GENO_TRBD2$FREQ_01 > 0.25 & Adaptive_GENO_TRBD2$FREQ_01 < 0.75]
D2_02_HOMO_FREQ <- Adaptive_GENO_TRBD2$FREQ_01[Adaptive_GENO_TRBD2$FREQ_01 < 0.25]

D2_01_HOMO_HETERO_BOUND <- (mean(D2_01_HOMO_FREQ)*sd(D2_HETERO_FREQ) + mean(D2_HETERO_FREQ)*sd(D2_01_HOMO_FREQ)) / (sd(D2_01_HOMO_FREQ)+sd(D2_HETERO_FREQ))
D2_02_HOMO_HETERO_BOUND <- (mean(D2_02_HOMO_FREQ)*sd(D2_HETERO_FREQ) + mean(D2_HETERO_FREQ)*sd(D2_02_HOMO_FREQ)) / (sd(D2_02_HOMO_FREQ)+sd(D2_HETERO_FREQ))

adapt_d2_error_graph <- ggplot(Adaptive_GENO_TRBD2, aes(x=FREQ_01)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.025), colour="black", fill="grey") +
  geom_vline(xintercept=c(D2_02_HOMO_HETERO_BOUND, D2_01_HOMO_HETERO_BOUND)) +
  xlab("TRBD2*01 fraction") + ylab("Number of individuals")+
  theme_classic() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=12), axis.text.x = element_text(size=12),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14))

adapt_d2_error_graph

require(grid)
title.grob <- textGrob(
  label = "C.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

adapt_d2_error_graph <- arrangeGrob(adapt_d2_error_graph, top = title.grob)


# ggsave(paste0(figure_folder, "DS4_TRBD2_allele_distribution.png"), adapt_d2_error_graph, height = 5)
# ggsave(paste0(figure_folder, "DS4_TRBD2_allele_distribution.pdf"), adapt_d2_error_graph, height = 5)

##################################################################################################################################
############################ DS3 (Cancer) TRBD2 ERROR GRAPH #########################################################################
##################################################################################################################################

names(DS3_genotypes) <- toupper(names(DS3_genotypes))
BIOMED2_GENO_TRBD2 <- DS3_genotypes[DS3_genotypes$GENE=="TRBD2",]
BIOMED2_GENO_TRBD2 <- BIOMED2_GENO_TRBD2[!grepl("[a-d]$", BIOMED2_GENO_TRBD2$SUBJECT),]
BIOMED2_GENO_TRBD2 <- BIOMED2_GENO_TRBD2[!grepl("hem", BIOMED2_GENO_TRBD2$SUBJECT),]

BIOMED2_GENO_TRBD2 <- separate(BIOMED2_GENO_TRBD2, "ALLELES", paste("ALLELES", 1:2, sep = "_"), sep = ",")
BIOMED2_GENO_TRBD2 <- separate(BIOMED2_GENO_TRBD2, "COUNTS", paste("COUNTS", 1:2, sep = "_"), sep = ",")

BIOMED2_GENO_TRBD2$COUNT_01 <- 0
BIOMED2_GENO_TRBD2$COUNT_01[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_1)] <- BIOMED2_GENO_TRBD2$COUNTS_1[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_1)] 
BIOMED2_GENO_TRBD2$COUNT_01[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_2)] <- BIOMED2_GENO_TRBD2$COUNTS_2[grepl("01", BIOMED2_GENO_TRBD2$ALLELES_2)] 

BIOMED2_GENO_TRBD2$COUNT_01 <- as.numeric(BIOMED2_GENO_TRBD2$COUNT_01)
BIOMED2_GENO_TRBD2$FREQ_01 <- BIOMED2_GENO_TRBD2$COUNT_01 / BIOMED2_GENO_TRBD2$TOTAL

summary(BIOMED2_GENO_TRBD2$TOTAL)

D2_HOMO02_BOUND <- (D2_02_HOMO_HETERO_BOUND + mean(D2_02_HOMO_FREQ)) / 2
D2_HETERO_1_BOUND <- (D2_02_HOMO_HETERO_BOUND + mean(D2_HETERO_FREQ)) / 2
D2_HETERO_2_BOUND <- (D2_01_HOMO_HETERO_BOUND + mean(D2_HETERO_FREQ)) / 2
D2_HOMO01_BOUND <- (D2_01_HOMO_HETERO_BOUND + mean(D2_01_HOMO_FREQ)) / 2

# d <- ggplot(ALL_GENO_TRBD2, aes(x=FREQ_01)) + 
#   geom_density(alpha=.2, fill="#FF6666") + ggtitle("TRBD2 Biomed2 dataset")

biomed_d2_error_graph <- ggplot(BIOMED2_GENO_TRBD2, aes(x=FREQ_01)) + 
  geom_rect(data=NULL,aes(xmin=D2_HOMO02_BOUND,xmax=D2_HETERO_1_BOUND,ymin=-Inf,ymax=Inf),
            fill="lightgrey")+
  geom_rect(data=NULL,aes(xmin=D2_HETERO_2_BOUND,xmax=D2_HOMO01_BOUND,ymin=-Inf,ymax=Inf),
            fill="lightgrey")+
  geom_histogram(breaks = seq(0, 1, by = 0.025), colour="black", fill="darkgrey") +
  geom_vline(xintercept=c(D2_02_HOMO_HETERO_BOUND, D2_01_HOMO_HETERO_BOUND)) +
  geom_vline(xintercept=c(D2_HOMO02_BOUND, D2_HETERO_1_BOUND, D2_HETERO_2_BOUND, D2_HOMO01_BOUND), linetype='dashed') +
  xlab("TRBD2*01 fraction") + ylab("Number of individuals")+
  theme_classic() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=12), axis.text.x = element_text(size=12),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14))


biomed_d2_error_graph

require(grid)
title.grob <- textGrob(
  label = "B.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

biomed_d2_error_graph <- arrangeGrob(biomed_d2_error_graph, top = title.grob)

# ggsave(paste0(figure_folder, "DS3_TRBD2_allele_distribution.png"), biomed_d2_error_graph, height = 5)
# ggsave(paste0(figure_folder, "DS3_TRBD2_allele_distribution.pdf"), biomed_d2_error_graph, height = 5)

##################################################################################################################################
############################ DS3 (Cancer) TRBJ1-6*01 frequency graph #############################################################
##################################################################################################################################

BIOMED2_GENO_TRBJ1_6 <- DS3_genotypes[DS3_genotypes$GENE=="TRBJ1-6",]
BIOMED2_GENO_TRBJ1_6 <- BIOMED2_GENO_TRBJ1_6[!grepl("[a-d]$", BIOMED2_GENO_TRBJ1_6$SUBJECT),]

BIOMED2_GENO_TRBJ1_6 <- separate(BIOMED2_GENO_TRBJ1_6, "ALLELES", paste("ALLELES", 1:2, sep = "_"), sep = ",")
BIOMED2_GENO_TRBJ1_6 <- separate(BIOMED2_GENO_TRBJ1_6, "COUNTS", paste("COUNTS", 1:2, sep = "_"), sep = ",")

BIOMED2_GENO_TRBJ1_6$COUNT_01 <- 0
BIOMED2_GENO_TRBJ1_6$COUNT_01[grepl("01", BIOMED2_GENO_TRBJ1_6$ALLELES_1)] <- BIOMED2_GENO_TRBJ1_6$COUNTS_1[grepl("01", BIOMED2_GENO_TRBJ1_6$ALLELES_1)] 
BIOMED2_GENO_TRBJ1_6$COUNT_01[grepl("01", BIOMED2_GENO_TRBJ1_6$ALLELES_2)] <- BIOMED2_GENO_TRBJ1_6$COUNTS_2[grepl("01", BIOMED2_GENO_TRBJ1_6$ALLELES_2)] 

BIOMED2_GENO_TRBJ1_6$COUNT_01 <- as.numeric(BIOMED2_GENO_TRBJ1_6$COUNT_01)
BIOMED2_GENO_TRBJ1_6$FREQ_01 <- BIOMED2_GENO_TRBJ1_6$COUNT_01 / BIOMED2_GENO_TRBJ1_6$TOTAL

summary(BIOMED2_GENO_TRBJ1_6$TOTAL)

biomed_j1_6_graph <- ggplot(BIOMED2_GENO_TRBJ1_6, aes(x=FREQ_01)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.025), colour="black", fill="darkgrey") +
  xlab("TRBJ1-6*01 fraction") + ylab("Number of individuals")+
  theme_classic()  + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=12), axis.text.x = element_text(size=12),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14))


biomed_j1_6_graph

require(grid)
title.grob <- textGrob(
  label = "D.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

biomed_j1_6_graph <- arrangeGrob(biomed_j1_6_graph, top = title.grob)

# ggsave(paste0(figure_folder, "DS3_TRBJ1_6_FREQ_01_02.pdf"), biomed_j1_6_graph, height = 5)
# ggsave(paste0(figure_folder, "DS3_TRBJ1_6_FREQ_01_02.png"), biomed_j1_6_graph, height = 5)


#####################################################################################################################
############################ Save figure ############################################################################
#####################################################################################################################


library(gridExtra)
# Merge all 4 graphs into 1 figure
g <- grid.arrange(
  hcv_d2_error_graph, biomed_d2_error_graph, adapt_d2_error_graph, biomed_j1_6_graph, 
  layout_matrix = rbind(c(1, 3),
                        c(2, 4))
)

ggsave(paste0(figure_folder, "D2_ERROR_AND_J1_6_RELATION.pdf"), g, width = 10, height = 8)
ggsave(paste0(figure_folder, "D2_ERROR_AND_J1_6_RELATION.png"), g, width = 10, height = 8)

