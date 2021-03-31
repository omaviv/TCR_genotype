library(ggplot2)
library(stringr)
library(dplyr)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_1/")

unique_seq_df <- read.delim(paste0(required_files_folder, "unique_sequences_per_ds.tsv"), sep = "\t", stringsAsFactors = F)

#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

unique_seq_plot <- ggplot(unique_seq_df, aes(x = dataset, y = log10(unique_reads))) +  
  geom_boxplot(aes(fill = factor(dataset))) +
  # geom_jitter(aes(colour = factor(dataset)), size = 1.5) +
  labs(fill="Data-set") +
  ylab("Log10(#unique sequences)") + xlab("Data-set")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16), axis.text = element_text(size = 16))

ggsave(paste0(figure_folder, "unique_sequences_boxplot.pdf"), unique_seq_plot)
# ggsave(paste0(figure_folder, "unique_sequences_boxplot.png"), unique_seq_plot)
