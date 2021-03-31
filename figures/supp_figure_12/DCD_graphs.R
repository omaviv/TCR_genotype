library(dplyr)
library(tidyr)
library(tidyverse)
library(seqinr)
library(gridExtra)
library(grid)

##################################################################################################################################################
############################ Load required files #################################################################################################
##################################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_12/")


genes_usage <- read.delim(paste0(required_files_folder, "HCV_DELETIONS.tab"), sep = "\t", stringsAsFactors = F)

##################################################################################################################################################
###########################  Figure  #############################################################################################################
##################################################################################################################################################

deleted_rows <- genes_usage[!is.na(genes_usage$col)& genes_usage$col=="Deletion",]
subjects_with_del <- genes_usage[!is.na(genes_usage$col)& genes_usage$col=="Deletion" & genes_usage$GENE=="TRBV4-3","SUBJECT"]


TRBV6_23_FREQ <- genes_usage[genes_usage$GENE=="TRBV6-23",]
TRBV6_23_FREQ$GROUP <- "Carry  of V4-3 and V3-2"
TRBV6_23_FREQ[TRBV6_23_FREQ$SUBJECT %in% subjects_with_del,"GROUP"] <- "V4-3 and V3-2 are deleted on"


p <- ggplot(TRBV6_23_FREQ, aes(x = GROUP, y = FREQ)) +  
  geom_boxplot(alpha = 0.80, fill=NA) +
  # geom_point(aes(colour = factor(HAPLO)), size = 1) + 
  geom_point(aes(colour = factor(GROUP)), size = 1) + 
  ggtitle("V6-2/V6-3") + xlab("Individual groups") + ylab("V6-2/V6-3 Usage") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none")


ggsave(paste0(figure_folder, "TRBV6-23_FREQ.pdf"), plot = p)
# ggsave(paste0(figure_folder, "TRBV6-23_FREQ.png"), plot = p)

