library(ggplot2)
library(gplots)
library(dplyr)
library(data.table)
library(fastmatch)
library(reshape2)
library(stringdist)
library(seqinr)
library(stringr)
library(stringi)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_8/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

upstream_v_seqs_df <- read.delim(paste0(required_files_folder, "5UTR_unfilter_seqs.tab"), header = T, sep = "\t", stringsAsFactors = F)
TRBV_LEADER <- read.fasta(paste0(required_files_folder, "TRBV_leader.fasta"),as.string = TRUE)

alleles2show <- c("TRBV20-1*01", "TRBV20-1*02", "TRBV4-3*01")

#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

upstream_v_seqs_df$ALLELE <- gsub("_[0-9]+[A-Z]+[0-9]+", "", upstream_v_seqs_df$ALLELE)
upstream_v_seqs_df <- upstream_v_seqs_df[(upstream_v_seqs_df$CLUSTER_COUNT > 2) | upstream_v_seqs_df$FRAC > 0.3,]

upstream_v_seqs_df <- upstream_v_seqs_df[upstream_v_seqs_df$ALLELE %in% alleles2show,]

alleles <- unique(upstream_v_seqs_df$ALLELE)

upstream_v_seqs_df <- upstream_v_seqs_df %>% select(CONS, GENE, ALLELE)
upstream_v_seqs_df$ALLELE <- paste0(upstream_v_seqs_df$ALLELE, "_consensus")

names(TRBV_LEADER) <- sapply(strsplit(names(TRBV_LEADER),"|",fixed = T),"[",2)
TRBV_LEADER <- toupper(unlist(TRBV_LEADER))
TRBV_LEADER <- TRBV_LEADER[names(TRBV_LEADER) %in% alleles2show]

imgt_alleles <- data.frame(CONS=stri_reverse(TRBV_LEADER), ALLELE=names(TRBV_LEADER), stringsAsFactors = F)
imgt_alleles$GENE <- sapply(strsplit(imgt_alleles$ALLELE, "*", fixed = T), "[", 1)
imgt_alleles$ALLELE <- paste0(imgt_alleles$ALLELE, "_imgt_ref")

upstream_v_seqs_df <- rbind(upstream_v_seqs_df, imgt_alleles)
## Filter 
# UTR.db.flt <- UTR.db[(UTR.db$FRAC>=0.1 & UTR.db$CLUSTER_COUNT>=4),]

UTR.db.flt <- upstream_v_seqs_df
###

UTR.db.flt$FAMILY <- sapply(strsplit(UTR.db.flt$GENE,'-',fixed=T),'[',1)

for.clust <- UTR.db.flt$CONS
names(for.clust) <- UTR.db.flt$ALLELE
dist.mat <- as.matrix(stringdistmatrix(for.clust,method = 'dl',useNames = F,nthread = 40))

clust.tmp <- hclust(dist(dist.mat))


################ PLOT 5UTR SEQUENCE


UTR.db.flt.rev <- UTR.db.flt


UTR.db.flt.mlt.rev <- do.call(rbind,lapply(1:nrow(UTR.db.flt.rev),function(i){
  x <- do.call(rbind,lapply(1:nchar(UTR.db.flt.rev$CONS[i]),function(j){UTR.db.flt.rev[i,]}));
  x$NT <- s2c(UTR.db.flt.rev$CONS[i]); 
  x$POS <- as.numeric(1:nchar(UTR.db.flt.rev$CONS[i]));
  return(x)}))

hold_UTR.db.flt.mlt.rev <- UTR.db.flt.mlt.rev
UTR.db.flt.mlt.rev <- hold_UTR.db.flt.mlt.rev

max(UTR.db.flt.mlt.rev$POS)
UTR.db.flt.mlt.rev$POS_REV <- max(UTR.db.flt.mlt.rev$POS)-UTR.db.flt.mlt.rev$POS+1

allele2gap <- "TRBV20-1*01_consensus"
position_2_open_gap <- 12
gap_len <- 33
UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS>=position_2_open_gap] <- UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS>=position_2_open_gap] + gap_len
UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS<position_2_open_gap] <- UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS<position_2_open_gap] + gap_len
temp_row <- UTR.db.flt.mlt.rev[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS==position_2_open_gap-1,]
temp_row$NT <- "gap"
for (pos in position_2_open_gap:(position_2_open_gap+gap_len-1)) {
  temp_row$POS <- pos
  UTR.db.flt.mlt.rev <- rbind(UTR.db.flt.mlt.rev, temp_row)
}

allele2gap <- "TRBV20-1*02_consensus"
position_2_open_gap <- 12
gap_len <- 3
UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS>=position_2_open_gap] <- UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS>=position_2_open_gap] + gap_len
UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS<position_2_open_gap] <- UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS<position_2_open_gap] + gap_len
temp_row <- UTR.db.flt.mlt.rev[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS==position_2_open_gap-1,]
temp_row$NT <- "gap"
for (pos in position_2_open_gap:(position_2_open_gap+gap_len-1)) {
  temp_row$POS <- pos
  UTR.db.flt.mlt.rev <- rbind(UTR.db.flt.mlt.rev, temp_row)
}


allele2gap <- "TRBV4-3*01_consensus"
position_2_open_gap <- 9
gap_len <- 9
UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS>=position_2_open_gap] <- UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS>=position_2_open_gap] + gap_len
UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS<position_2_open_gap] <- UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS<position_2_open_gap] + gap_len
temp_row <- UTR.db.flt.mlt.rev[UTR.db.flt.mlt.rev$ALLELE == allele2gap & UTR.db.flt.mlt.rev$POS==position_2_open_gap-1,]
temp_row$NT <- "gap"
for (pos in position_2_open_gap:(position_2_open_gap+gap_len-1)) {
  temp_row$POS <- pos
  UTR.db.flt.mlt.rev <- rbind(UTR.db.flt.mlt.rev, temp_row)
}



UTR.db.flt.mlt.rev$NT <- factor(UTR.db.flt.mlt.rev$NT,levels=c('A','C','G','T','N','gap'))

# UTR.db.flt.mlt.rev$ALLELE <- gsub("TRB", "", UTR.db.flt.mlt.rev$ALLELE)
UTR.db.flt.mlt.rev$ALLELE <- factor(UTR.db.flt.mlt.rev$ALLELE, levels=sort(unique(UTR.db.flt.mlt.rev$ALLELE), decreasing = TRUE))

FAM.COL <- c('TRBV1'='red','TRBV2'='blue','TRBV3'='darkgreen','TRBV4'='purple','TRBV5'='orange','TRBV6'='cyan','TRBV7'='black')
color.panel2 <- c('A'='#e69f00','C'='#56b4e9','G'='#009e73','T'='#f0e442', 'N'='grey', 'gap'='black')


pdf(paste0(figure_folder, "5utr_against_ref.pdf"),height = 8,width = 10)

break_lim <- round(max(UTR.db.flt.mlt.rev$POS)/10)*10
ggplot(UTR.db.flt.mlt.rev, aes(x=(-POS), y=ALLELE, fill=NT)) + 
  geom_tile(colour="white") +
  scale_fill_manual(values = color.panel2) +  
  theme(axis.title.y=element_blank(),axis.text.x=element_text(size=12),
        axis.text.y = element_text(size=8,color = "black"),
        legend.title=element_blank(),legend.position="left",
        axis.title.x = element_text(size = 14))+
  scale_y_discrete(position = "right") + 
  scale_x_continuous(breaks=seq(0,-break_lim,-10), expand = c(0, 0))+
  xlab('Position')


dev.off()

