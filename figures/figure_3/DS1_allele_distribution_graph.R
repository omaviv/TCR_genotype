library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(wesanderson)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_3/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")
sources_folder <- paste0(project_folder, "figures/sources/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

# load sources
source(paste0(sources_folder, "genotype_heatmap/internal_functions.R"))
source(paste0(sources_folder, "genotype_heatmap/genoHeatmap.R"))

# load DS1 genotypes after filtering unreliable unknown alleles
genotypes <- read.delim(paste0(required_files_folder, "HCV_Genotypes_filtered.tab"), sep = "\t", stringsAsFactors = F)

###################################################################################################################################
############################ Genrate allele distribution graphs ###################################################################
###################################################################################################################################

genotypes <- genotypes[!genotypes$GENE %in% PSEUDO[["TRB"]],]
genotypes$GENOTYPED_ALLELES <- gsub("_[0-9]+[A-Z]+[0-9]+", "", genotypes$GENOTYPED_ALLELES)

homozygous <- genotypes[!grepl(",", genotypes$GENOTYPED_ALLELES),]
homozygous <- homozygous %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N = n()*2)

heterozygous <- genotypes[grepl(",", genotypes$GENOTYPED_ALLELES),]
heterozygous <- heterozygous %>% separate_rows(GENOTYPED_ALLELES, sep = ",")
heterozygous <- heterozygous %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N = n())

allele_dis <- rbind(homozygous, heterozygous)
# unique(gsub("_[0-9]+$", "", allele_dis$GENOTYPED_ALLELES)) <- "03_313GCCATCAGTGAGTC326"

# allele_dis$GENOTYPED_ALLELES[grepl("_[0-9]", allele_dis$GENOTYPED_ALLELES)] <- sapply(
#   strsplit(allele_dis$GENOTYPED_ALLELES[grepl("_[0-9]", allele_dis$GENOTYPED_ALLELES)], "_", fixed = T), "[", 1)
allele_dis <- allele_dis %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N = sum(N))
allele_dis$TOTAL_GENE <- unlist(lapply(allele_dis$GENE, function(gene) {sum(allele_dis$N[allele_dis$GENE == gene])}))
allele_dis$FREQ <- allele_dis$N / allele_dis$TOTAL_GENE


allele_dis$GENOTYPED_ALLELES <- gsub("Deletion","Del",allele_dis$GENOTYPED_ALLELES)
allele_palette <- allelePalette(allele_dis$GENOTYPED_ALLELES)


novel <- data.frame(Novel=grep('[_]',allele_dis$GENOTYPED_ALLELES,value=T),
                    Base = sapply(grep('[_]',allele_dis$GENOTYPED_ALLELES,value=T),
                                  function(x) strsplit(x,'[_]')[[1]][1])) %>%
  distinct() %>% group_by(.data$Base) %>% dplyr::mutate(n = dplyr::n())

dagger = "^"
id <- grep('^[0-9]+[_].',names(allele_palette$transper))
allele_palette$transper[id] <- 1
new_allele <- paste0(dagger,1:length(id),'-',allele_palette$AlleleCol[id])
names(new_allele) <-allele_palette$AlleleCol[id]

for(i in 1:length(new_allele)){
  allele <- names(new_allele)[i]
  allele_dis$GENOTYPED_ALLELES[grep(allele,allele_dis$GENOTYPED_ALLELES)] <- new_allele[i]
}
allele_palette$AlleleCol[id] <- new_allele
names(allele_palette$transper)[id] <- new_allele
novel <- T

allele_dis$GENOTYPED_ALLELES <- factor(allele_dis$GENOTYPED_ALLELES, levels = allele_palette$AlleleCol)
allele_dis$TEXT <- sapply(strsplit(as.character(allele_dis$GENOTYPED_ALLELES), "-", fixed = T), "[", 1)
allele_dis$TEXT[grepl("^0[0-9]|D", allele_dis$TEXT)] <- ""


allele_dis$pos <- 0
for (gene in allele_dis$GENE) {
  temp <- allele_dis[allele_dis$GENE == gene,]
  total_freq <- 0
  for (i in 1:nrow(temp)) {
    temp$pos[[i]] <- total_freq + (temp$FREQ[[i]] / 2)
    total_freq <- total_freq + temp$FREQ[[i]]
  }
  allele_dis$pos[allele_dis$GENE == gene] <- temp$pos
}

allele_dis$pos <- 1 - allele_dis$pos


v_allele_dis <- allele_dis[grepl("V", allele_dis$GENE),]
v_allele_dis <- v_allele_dis[!grepl("OR", v_allele_dis$GENE),]
v_allele_dis$GENE <- gsub("TRB", "", v_allele_dis$GENE)


TRBV_LOC <- GENE.loc[["TRB"]][grepl("V", GENE.loc[["TRB"]])]
TRBV_LOC <- TRBV_LOC[!TRBV_LOC %in% PSEUDO[["TRB"]]]

TRBV_LOC <- gsub("TRB", "",TRBV_LOC)
TRBV_LOC <- TRBV_LOC[TRBV_LOC %in% unique(v_allele_dis$GENE)]


v_allele_dis_graph <- ggplot(v_allele_dis, aes(x=GENE, y=FREQ, fill=GENOTYPED_ALLELES)) + 
  geom_bar(stat="identity", colour="white") + 
  geom_text(aes_string(label = "TEXT", x = "GENE", y = "pos",text = "TEXT"), size=3) +
  scale_y_continuous(trans = "reverse", expand = c(0, 0)) + 
  scale_x_discrete(position = "top", limits=TRBV_LOC) +
  scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE) +
  labs(fill="Alleles") + ylab("Distribution") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14), axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "none", 
        axis.title.y = element_text(size = 26))

# D and J allele distribution
dj_allele_dis <- allele_dis[!grepl("V", allele_dis$GENE),]
dj_allele_dis$GENE <- gsub("TRB", "", dj_allele_dis$GENE)

TRBDJ_LOC <- GENE.loc[["TRB"]][!grepl("V", GENE.loc[["TRB"]])]
TRBDJ_LOC <- TRBDJ_LOC[!TRBDJ_LOC %in% PSEUDO[["TRB"]]]

# TRBDJ_LOC <- TRBDJ_LOC[!grepl("P", TRBDJ_LOC)]

dj_allele_dis_graph <- ggplot(dj_allele_dis, aes(x=GENE, y=FREQ, fill=GENOTYPED_ALLELES)) + 
  geom_bar(stat="identity", colour="white") + 
  geom_text(aes_string(label = "TEXT", x = "GENE", y = "pos",text = "TEXT"), size=3) +
  scale_y_continuous(trans = "reverse", expand = c(0, 0)) + 
  scale_x_discrete(position = "top", limits=gsub("TRB", "",TRBDJ_LOC)) +
  scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE) +
  labs(fill="Alleles") + ylab("Distribution") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14), axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "none", 
        axis.title.y = element_text(size = 26))


#####################################################################################################################
############################ Save graphs ############################################################################
#####################################################################################################################

ggsave(paste0(figure_folder, "DS1_v_alleles_distribution.pdf"), v_allele_dis_graph, height = 6, width = 18.5)
# ggsave(paste0(figure_folder, "DS1_v_alleles_distribution.png"), v_allele_dis_graph, height = 6, width = 18.5)

ggsave(paste0(figure_folder, "DS1_dj_alleles_distribution.pdf"), dj_allele_dis_graph, height = 6, width = 5.6)
# ggsave(paste0(figure_folder, "DS1_dj_alleles_distribution.png"), dj_allele_dis_graph, height = 6, width = 5)


