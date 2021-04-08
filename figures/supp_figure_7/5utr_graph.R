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
figure_folder <- paste0(project_folder, "figures/supp_figure_7/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

upstream_v_seqs_df <- read.delim(paste0(required_files_folder, "5UTR_unfilter_seqs.tab"), header = T, sep = "\t", stringsAsFactors = F)
TRBV_LEADER <- read.fasta(paste0(required_files_folder, "TRBV_leader.fasta"),as.string = TRUE)

# Change the names of the collapsed genes
# TRBV6-23 > TRBV6-2/TRBV6-3
upstream_v_seqs_df$GENE <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", upstream_v_seqs_df$GENE)
upstream_v_seqs_df$ALLELE <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", upstream_v_seqs_df$ALLELE)
names(TRBV_LEADER) <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", names(TRBV_LEADER))

#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

upstream_v_seqs_df$ALLELE <- gsub("_[0-9]+[A-Z]+[0-9]+", "", upstream_v_seqs_df$ALLELE)
upstream_v_seqs_df <- upstream_v_seqs_df[(upstream_v_seqs_df$CLUSTER_COUNT > 2) | upstream_v_seqs_df$FRAC > 0.3,]

alleles <- unique(upstream_v_seqs_df$ALLELE)

for (allele in alleles) {
  temp <- upstream_v_seqs_df[upstream_v_seqs_df$ALLELE==allele,]
  n_cons <- nrow(temp)
  if (nrow(temp) > 1) {
    for (i in 1:10){
      cut_seq <- substr(temp$CONS, 1, (str_length(temp$CONS)[[1]] - i))
      if (length(unique(cut_seq)) < n_cons) {
        upstream_v_seqs_df$CONS[upstream_v_seqs_df$ALLELE == allele] <- sapply(upstream_v_seqs_df$CONS[upstream_v_seqs_df$ALLELE == allele], function(seq) {
          substr(seq, 1, str_length(seq) - i)
        })
        n_cons <- length(unique(cut_seq))[[1]]
      }
    }
  }
}


upstream_v_seqs_df <- upstream_v_seqs_df %>% select(CONS, GENE, SAMP, CLUSTER_COUNT, ALLELE, FRAC)

upstream_v_seqs_df <- upstream_v_seqs_df %>% group_by(CONS, GENE, ALLELE) %>% dplyr::summarise(SAMP = paste(SAMP, collapse = ","), FRAC = sum(FRAC))
upstream_v_seqs_df$SAMP <- sapply(upstream_v_seqs_df$SAMP, function(x) {paste(unique(strsplit(x, ",", fixed = T)[[1]]), collapse=",")})
upstream_v_seqs_df$CLUSTER_COUNT <- str_count(upstream_v_seqs_df$SAMP, ",") + 1

## Filter 
# UTR.db.flt <- UTR.db[(UTR.db$FRAC>=0.1 & UTR.db$CLUSTER_COUNT>=4),]

UTR.db.flt <- upstream_v_seqs_df
###

UTR.db.flt$FAMILY <- sapply(strsplit(UTR.db.flt$GENE,'-',fixed=T),'[',1)

names(TRBV_LEADER) <- sapply(strsplit(names(TRBV_LEADER),"|",fixed = T),"[",2)
TRBV_LEADER <- toupper(unlist(TRBV_LEADER))

UTR.db.flt$as_germ <- unlist(lapply(1:nrow(UTR.db.flt), function(i){
  allele <- sapply(strsplit(UTR.db.flt$ALLELE[[i]],"_",fixed = T),"[",1);
  if (!allele %in% names(TRBV_LEADER)) {return(FALSE)}
  allele_seq <- stri_reverse(TRBV_LEADER[[allele]]);
  return(grepl(allele_seq, UTR.db.flt$CONS[[i]], fixed = T))
}))

UTR.db.flt <- UTR.db.flt[order(UTR.db.flt$ALLELE, -UTR.db.flt$as_germ),]


if(length(unique(UTR.db.flt$ALLELE)) != length(UTR.db.flt$ALLELE)){
  dup.samp <- UTR.db.flt$ALLELE[duplicated(UTR.db.flt$ALLELE)]
  for(dup in dup.samp){
    UTR.db.flt$ALLELE[UTR.db.flt$ALLELE==dup] <- paste0(UTR.db.flt$ALLELE[UTR.db.flt$ALLELE==dup],'_',
                                                        1:length(UTR.db.flt$ALLELE[UTR.db.flt$ALLELE==dup]))
  }
}

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

UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == "TRBV23-1*01_1" & UTR.db.flt.mlt.rev$POS>=89] <- UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == "TRBV23-1*01_1" & UTR.db.flt.mlt.rev$POS>=89] + 5
UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == "TRBV23-1*01_1" & UTR.db.flt.mlt.rev$POS<89] <- UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == "TRBV23-1*01_1" & UTR.db.flt.mlt.rev$POS<89] + 5
temp_row <- UTR.db.flt.mlt.rev[UTR.db.flt.mlt.rev$ALLELE == "TRBV23-1*01_1" & UTR.db.flt.mlt.rev$POS==88,]
temp_row$NT <- "gap"
for (pos in 89:93) {
  temp_row$POS <- pos
  UTR.db.flt.mlt.rev <- rbind(UTR.db.flt.mlt.rev, temp_row)
}

UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == "TRBV7-7*01_C315T_1" & UTR.db.flt.mlt.rev$POS>=9] <- UTR.db.flt.mlt.rev$POS[UTR.db.flt.mlt.rev$ALLELE == "TRBV7-7*01_C315T_1" & UTR.db.flt.mlt.rev$POS>=9] + 15
UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == "TRBV7-7*01_C315T_1" & UTR.db.flt.mlt.rev$POS<9] <- UTR.db.flt.mlt.rev$POS_REV[UTR.db.flt.mlt.rev$ALLELE == "TRBV7-7*01_C315T_1" & UTR.db.flt.mlt.rev$POS<9] + 15
temp_row <- UTR.db.flt.mlt.rev[UTR.db.flt.mlt.rev$ALLELE == "TRBV7-7*01_C315T_1" & UTR.db.flt.mlt.rev$POS==8,]
temp_row$NT <- "gap"
for (pos in 9:23) {
  temp_row$POS <- pos
  UTR.db.flt.mlt.rev <- rbind(UTR.db.flt.mlt.rev, temp_row)
}

UTR.db.flt.mlt.rev$NT <- factor(UTR.db.flt.mlt.rev$NT,levels=c('A','C','G','T','N','gap'))

# UTR.db.flt.mlt.rev$ALLELE <- gsub("TRB", "", UTR.db.flt.mlt.rev$ALLELE)
UTR.db.flt.mlt.rev$ALLELE <- factor(UTR.db.flt.mlt.rev$ALLELE, levels=sort(unique(UTR.db.flt.mlt.rev$ALLELE), decreasing = TRUE))

FAM.COL <- c('TRBV1'='red','TRBV2'='blue','TRBV3'='darkgreen','TRBV4'='purple','TRBV5'='orange','TRBV6'='cyan','TRBV7'='black')
color.panel2 <- c('A'='#e69f00','C'='#56b4e9','G'='#009e73','T'='#f0e442', 'N'='grey', 'gap'='black')


pdf(paste0(figure_folder, "5utr_db_geno.pdf"),height = 18,width = 14)

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

