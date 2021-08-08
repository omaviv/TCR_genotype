library(dplyr)
library(tidyr)
library(ggplot2)
library(seqinr)
library(stringr)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_2/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

Adaptive_V_start <- read.delim(paste0(required_files_folder, "Adaptive_V_start.tab"), header=T, sep = "\t", stringsAsFactors = F)
TRBV_GERM <- read.fasta(paste0(references_folder, "Adaptive_TRBV.fasta"),as.string = TRUE)

#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

TRBV_GERM <- toupper(unlist(TRBV_GERM))

allele_ndots <- list()
for (allele in names(TRBV_GERM)) {
  allele_ndots[allele] <- str_count(TRBV_GERM[[allele]], pattern = "[.]")
}

Adaptive_V_start_novel <- Adaptive_V_start[grepl("_", Adaptive_V_start$v_call, fixed = T),]

Adaptive_V_start_novel$v_germline_ndots <- unlist(lapply(Adaptive_V_start_novel$v_call, function(v) {
  paste(unique(allele_ndots[names(allele_ndots) %in% strsplit(strsplit(v, ",", fixed = T)[[1]], "_", fixed = T)[[1]]]), collapse = ",")
}))

Adaptive_V_start_novel$v_germline_ndots <- as.integer(Adaptive_V_start_novel$v_germline_ndots)
Adaptive_V_start_novel$v_ref_start <- Adaptive_V_start_novel$germline_ndots - Adaptive_V_start_novel$v_germline_ndots + 1

Adaptive_V_start <- rbind(Adaptive_V_start[!grepl("_", Adaptive_V_start$v_call, fixed = T),], Adaptive_V_start_novel)
Adaptive_V_start <- Adaptive_V_start[complete.cases(Adaptive_V_start),]
Adaptive_V_start$start_pos <- Adaptive_V_start$v_ref_start - Adaptive_V_start$v_sequence_start

Adaptive_V_start$v_gene <- sapply(strsplit(Adaptive_V_start$v_call, "*", fixed = T), "[", 1)

v_germ_start_pos <- Adaptive_V_start %>% group_by(v_ref_start) %>% dplyr::summarise(N=sum(N))
v_germ_start_pos$Freq <- v_germ_start_pos$N / sum(v_germ_start_pos$N)

v_germ_start_pos_graph <- ggplot(v_germ_start_pos[v_germ_start_pos$v_ref_start <= 10,], aes(x=as.factor(v_ref_start-43), y=Freq))+
  geom_bar(stat="identity", colour="white") + 
  ylab("Frequency") + xlab("First covered position from the TRBV gene reference") +
  theme_classic()
v_germ_start_pos_graph

ggsave(paste0(figure_folder, "Adaptive_V_germ_start_pos_fig.pdf"), v_germ_start_pos_graph, width = 10)
# ggsave(paste0(figure_folder, "Adaptive_V_germ_start_pos_fig.png"), v_germ_start_pos_graph, width = 10)

