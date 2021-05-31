library(ggplot2)
library(gplots)
library(dplyr)
library(data.table)
library(fastmatch)
library(reshape2)
library(ggsignif)
library(gridExtra)


##################################################################################################################################
############################ Load required files #################################################################################
##################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_10/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

pairing_df <- read.delim(paste0(required_files_folder, "DS4_all_dj_pairs.tab"), sep = "\t", stringsAsFactors = F)
genotypes <- read.delim(paste0(required_files_folder, "Adaptive_All_Genotypes.tab"), sep = "\t", stringsAsFactors = F)

##################################################################################################################################################
################### Initial parameters and functions #############################################################################################
##################################################################################################################################################

require(grid)
set_title <- function(panel_symbol) {
  textGrob(
    label = panel_symbol,
    x = unit(0, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0)
}
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


TRBD2_geno <- genotypes[genotypes$GENE=="TRBD2",]
TRBD2_geno$GENOTYPED_ALLELES[grepl(",", TRBD2_geno$GENOTYPED_ALLELES)] <- "01,02"
TRBD2_geno$d2_geno <- TRBD2_geno$GENOTYPED_ALLELES
TRBD2_geno <- TRBD2_geno %>% select(SUBJECT, d2_geno)
names(TRBD2_geno) <- tolower(names(TRBD2_geno))
pairing_df <- merge(pairing_df, TRBD2_geno)

##################################################################################################################################################
################### TRBD gene multiple according to length #######################################################################################
##################################################################################################################################################

d_usage_by_len <- pairing_df %>% group_by(subject, d_gene, d_len) %>% dplyr::summarise(N=sum(number_of_comb))
total_d_len <- d_usage_by_len %>% group_by(subject, d_len) %>% dplyr::summarise(total=sum(N))
d_usage_by_len <- merge(d_usage_by_len, total_d_len, by=c("subject", "d_len"))
d_usage_by_len$fraction <- d_usage_by_len$N/d_usage_by_len$total
d_usage_by_len$d_gene <- gsub("TRB", "", d_usage_by_len$d_gene)

d_usage_by_len_graph <- ggplot(d_usage_by_len[d_usage_by_len$d_len <= 12,], aes(x=as.factor(d_len), y=fraction)) + 
  geom_boxplot(aes(fill = factor(d_gene))) +
  xlab("TRBD length") + ylab("Fraction") + labs(fill="Assignments") + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# ggsave(paste0(figure_folder, "DS4_TRBD_usage_by_length.pdf"), d_usage_by_len_graph)
# ggsave(paste0(figure_folder, "DS4_TRBD_usage_by_length.png"), d_usage_by_len_graph)

# ##################################################################################################################################################
# ################### TRBD gene multiple according to length D2*01 ################################################################################
# ##################################################################################################################################################
# d_usage_by_len <- pairing_df[pairing_df$d2_geno == "01",]
# d_usage_by_len <- d_usage_by_len %>% group_by(subject, d_gene, d_len) %>% dplyr::summarise(N=sum(number_of_comb))
# total_d_len <- d_usage_by_len %>% group_by(subject, d_len) %>% dplyr::summarise(total=sum(N))
# d_usage_by_len <- merge(d_usage_by_len, total_d_len, by=c("subject", "d_len"))
# d_usage_by_len$fraction <- d_usage_by_len$N/d_usage_by_len$total
# 
# d_usage_by_len_01_graph <- ggplot(d_usage_by_len[d_usage_by_len$d_len <= 12,], aes(x=as.factor(d_len), y=fraction)) + 
#   geom_boxplot(aes(fill = factor(d_gene)))+
#   xlab("D length") + ylab("Fraction") + labs(fill="Assignments") + 
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# # ggsave(paste0(figure_folder, "DS4_TRBD_usage_by_length_01.pdf"), d_usage_by_len_01_graph)
# # ggsave(paste0(figure_folder, "DS4_TRBD_usage_by_length_01.png"), d_usage_by_len_01_graph)
# 
# ##################################################################################################################################################
# ################### TRBD gene multiple according to length #######################################################################################
# ##################################################################################################################################################
# d_usage_by_len <- pairing_df[pairing_df$d2_geno == "02",]
# d_usage_by_len <- d_usage_by_len %>% group_by(subject, d_gene, d_len) %>% dplyr::summarise(N=sum(number_of_comb))
# total_d_len <- d_usage_by_len %>% group_by(subject, d_len) %>% dplyr::summarise(total=sum(N))
# d_usage_by_len <- merge(d_usage_by_len, total_d_len, by=c("subject", "d_len"))
# d_usage_by_len$fraction <- d_usage_by_len$N/d_usage_by_len$total
# 
# d_usage_by_len_02_graph <- ggplot(d_usage_by_len[d_usage_by_len$d_len <= 12,], aes(x=as.factor(d_len), y=fraction)) + 
#   geom_boxplot(aes(fill = factor(d_gene)))+
#   xlab("D length") + ylab("Fraction") + labs(fill="Assignments") + 
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# # ggsave(paste0(figure_folder, "DS4_TRBD_usage_by_length_02.pdf"), d_usage_by_len_02_graph)
# # ggsave(paste0(figure_folder, "DS4_TRBD_usage_by_length_02.png"), d_usage_by_len_02_graph)

##################################################################################################################################################
################### TRBJ1-TRBD2 pairing fraction according to length #############################################################################
##################################################################################################################################################

d2_j1_pairing <- pairing_df[grepl("J1", pairing_df$j_gene) & pairing_df$d_gene=="TRBD2",] 
d2_j1_pairing <- d2_j1_pairing %>% group_by(subject, d_gene, d_len, d2_geno) %>% dplyr::summarise(pair_frac=sum(pair_frac))

d2_j1_pairing_graph <- ggplot(d2_j1_pairing, aes(x=as.factor(d_len), y=pair_frac)) + 
  geom_boxplot(aes(fill=d2_geno))+
  xlab("TRBD length") + ylab("P(JTRB1|TRBD2)") + 
  guides(fill=guide_legend("TRBD2 genotype")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# ggsave(paste0(figure_folder, "DS4_TRBD2_J1_pairing_by_length.pdf"), d2_j1_pairing_graph)
# ggsave(paste0(figure_folder, "DS4_TRBD2_J1_pairing_by_length.png"), d2_j1_pairing_graph)

##################################################################################################################################################
################### Collapse the 4 plots into one figure #########################################################################################
##################################################################################################################################################

d_usage_by_len_graph <- arrangeGrob(d_usage_by_len_graph, top = set_title("a."))
# d_usage_by_len_01_graph <- arrangeGrob(d_usage_by_len_01_graph, top = set_title("b."))
# d_usage_by_len_02_graph <- arrangeGrob(d_usage_by_len_02_graph, top = set_title("c."))
d2_j1_pairing_graph <- arrangeGrob(d2_j1_pairing_graph, top = set_title("b."))


g <- grid.arrange(
  grobs = list(d_usage_by_len_graph,  d2_j1_pairing_graph),
  nrow = 1,
  layout_matrix = rbind(
    c(1, 2))
)

ggsave(paste0(figure_folder, "D_misidentification_by_length.pdf"), g, height = 4, width = 10)
# ggsave(paste0(figure_folder, "D_misidentification_by_length.png"), g, height = 4, width = 10)

