library(dplyr)
library(ggplot2)
library(stringr)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_3/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

TRBD2_GENO <- read.delim(paste0(required_files_folder, "HCV_Genotypes_const_2.tab"), header = T, sep = "\t", stringsAsFactors = F)
v_frequencies <- read.csv(paste0(required_files_folder, "HCV_DELETIONS.tab"), sep = "\t", stringsAsFactors = F)

# Change the names of the collapsed genes
# TRBV6-23 > TRBV6-2/TRBV6-3
GENE.loc[["TRB"]] <- gsub("TRBV6-23", "TRBV6-2/V6-3", GENE.loc[["TRB"]])
v_frequencies$GENE <- gsub("TRBV6-23", "TRBV6-2/V6-3", v_frequencies$GENE)

######################################## TRBD2 genotype 
TRBD2_GENO <- TRBD2_GENO[TRBD2_GENO$GENE == "TRBD2",]

d2_01_subjects <- TRBD2_GENO$SUBJECT[TRBD2_GENO$GENOTYPED_ALLELES == "01"]
d2_02_subjects <- TRBD2_GENO$SUBJECT[TRBD2_GENO$GENOTYPED_ALLELES == "02"]

####################################################################################################################################
##################################################### TRBV gene usage #############################################################
####################################################################################################################################


v_frequencies$FAMILY <- sapply(strsplit(v_frequencies$GENE, "-", fixed = T), "[", 1)

families_with_mul_genes <- c("TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7", "TRBV10", "TRBV11", "TRBV12")
v_frequencies$FAMILY[!v_frequencies$FAMILY %in% families_with_mul_genes] <- "Rest"

v_frequencies <- v_frequencies[!v_frequencies$GENE %in% PSEUDO[["TRB"]],]
TRBV_LOC <- TRBV_LOC[!TRBV_LOC %in% PSEUDO[["TRB"]]]

v_frequencies$GENE <- gsub("TRB", "", v_frequencies$GENE)
TRBV_LOC <- GENE.loc[["TRB"]][grepl("V", GENE.loc[["TRB"]])]
TRBV_LOC <- gsub("TRB", "",TRBV_LOC)
TRBV_LOC <- TRBV_LOC[TRBV_LOC %in% unique(v_frequencies$GENE)]
# TRBV_LOC <- TRBV_LOC[!grepl("V26", TRBV_LOC)]


v_frequencies$D2_GENO <- "01_02"
v_frequencies$D2_GENO[v_frequencies$SUBJECT %in% d2_01_subjects] <- "01"
v_frequencies$D2_GENO[v_frequencies$SUBJECT %in% d2_02_subjects] <- "02"

v_usage_graph <- ggplot(v_frequencies, aes(x = GENE, y = FREQ)) +  
  geom_boxplot(alpha = 0.80, fill=NA) +
  geom_jitter(aes(colour = factor(FAMILY)), size = 1.5) +
  # geom_jitter(aes(colour = factor(D2_GENO)), size = 1.5) +
  scale_x_discrete(limits=gsub("TRB", "", TRBV_LOC)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(color="Gene family") +
  ylab("Usage") +
  guides(colour=guide_legend(nrow = 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=28), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 24), legend.position=c(0.5,0.9))


####################################################################################################################################
##################################################### TRBDJ gene usage #############################################################
####################################################################################################################################


ALL_GENO <- read.delim("/home/aviv/master_project/genotypes/HCV_Genotypes_const_2.tab", header = T, sep = "\t", stringsAsFactors = F)

ALL_J_GENO <- ALL_GENO[grepl("TRBJ", ALL_GENO$GENE),]


total_trbj <- aggregate(ALL_J_GENO$TOTAL, by=list(SUBJECT=ALL_J_GENO$SUBJECT), FUN=sum)
names(total_trbj)[names(total_trbj) == "x"] <- "TOTAL_SEQS"
ALL_J_GENO <- merge(ALL_J_GENO, total_trbj, by = "SUBJECT")

ALL_J_GENO$FREQ <- ALL_J_GENO$TOTAL / ALL_J_GENO$TOTAL_SEQS


ALL_D_GENO <- ALL_GENO[grepl("TRBD", ALL_GENO$GENE),]

total_trbd <- aggregate(ALL_D_GENO$TOTAL, by=list(SUBJECT=ALL_D_GENO$SUBJECT), FUN=sum)
names(total_trbd)[names(total_trbd) == "x"] <- "TOTAL_SEQS"
ALL_D_GENO <- merge(ALL_D_GENO, total_trbd, by = "SUBJECT")

ALL_D_GENO$FREQ <- ALL_D_GENO$TOTAL / ALL_D_GENO$TOTAL_SEQS

ALL_DJ_GENO <- rbind(ALL_D_GENO, ALL_J_GENO)
ALL_DJ_GENO$FAMILY <- "1"
ALL_DJ_GENO$FAMILY[grepl("D2|J2", ALL_DJ_GENO$GENE)] <- "2"
TRBDJ_LOCUS <- c("TRBD1", "TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6",
                 "TRBD2", "TRBJ2-1", "TRBJ2-2", "TRBJ2-2P", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-6", "TRBJ2-7")

ALL_DJ_GENO <- ALL_DJ_GENO[!ALL_DJ_GENO$GENE %in% PSEUDO[["TRB"]],]
TRBDJ_LOCUS <- TRBDJ_LOCUS[!TRBDJ_LOCUS %in% PSEUDO[["TRB"]]]

ALL_DJ_GENO$GENE <- gsub("TRB", "", ALL_DJ_GENO$GENE)
TRBDJ_LOCUS <- gsub("TRB", "", TRBDJ_LOCUS)

ALL_DJ_GENO$D2_GENO <- "01_02"
ALL_DJ_GENO$D2_GENO[ALL_DJ_GENO$SUBJECT %in% d2_01_subjects] <- "01"
ALL_DJ_GENO$D2_GENO[ALL_DJ_GENO$SUBJECT %in% d2_02_subjects] <- "02"


dj_usage_graph <- ggplot(ALL_DJ_GENO, aes(x = GENE, y = FREQ)) +  
  geom_boxplot(alpha = 0.80, fill=NA) +
  # geom_jitter(aes(colour = factor(FAMILY)), size = 1.5)  + scale_color_hue(l=40, c=35) +
  geom_jitter(aes(colour = factor(D2_GENO)), size = 1.75)  + scale_color_hue(l=40, c=35) +
  # scale_color_hue(l=40, c=35) +
  scale_x_discrete(limits=TRBDJ_LOCUS) + scale_y_continuous(expand = c(0, 0)) +
  # labs(color="Gene group") +
  labs(color="TRBD2\ngenotype") +
  ylab("Usage") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=28), axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 28), legend.position=c(0.8,0.85))


#####################################################################################################################
############################ Save graphs ############################################################################
#####################################################################################################################

ggsave(paste0(figure_folder, "DS1_TRBDJ_usage_colored_by_D2_GENO.pdf"), dj_usage_graph, height = 8, width = 5.35)
# ggsave(paste0(figure_folder, "DS1_TRBDJ_usage_colored_by_D2_GENO.png"), dj_usage_graph, height = 8, width = 7)

ggsave(paste0(figure_folder, "DS1_TRBV_usage.pdf"), v_usage_graph, height = 8, width = 15)
# ggsave(paste0(figure_folder, "DS1_TRBV_usage.png"), v_usage_graph, height = 8, width = 15)
