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
library(tidyr)
library(gridExtra)

##################################################################################################################################
############################ Load required files #################################################################################
##################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/supp_figure_17/")

DATA <- read.delim(paste0(required_files_folder, "TRBV6_23_merged_makedb.tab"), header = T, sep = "\t", stringsAsFactors = F)
upstream_v_seqs_df <- read.delim(paste0(required_files_folder, "5UTR_unfilter_seqs.tab"), header = T, sep = "\t", stringsAsFactors = F)
TRBV_LEADER <- read.fasta(paste0(required_files_folder, "TRBV_leader.fasta"), as.string = TRUE)
HCV_ALL_GENO <- read.delim(paste0(required_files_folder, "HCV_Genotypes_const_2.tab"), header = T, sep = "\t", stringsAsFactors = F)


##################################################################################################################################
############################ Calculate the fraction of TRBV6-23 5' UTR over each subject #########################################
##################################################################################################################################

#Extract the 5 UTRs from sequences
DATA <- DATA[DATA$V_GERM_START_VDJ==1 & DATA$CONSCOUNT >= 2,]
DATA$UTR <- sapply(1:nrow(DATA),function(i){substr(start = 1,stop=DATA$V_SEQ_START[i]-1,DATA$SEQUENCE_INPUT[i])})
DATA$SEQ_NAME <- paste0(DATA$SUBJECT, "_", DATA$SEQUENCE_ID, "_", DATA$V_CALL)
DATA <- DATA[str_length(DATA$UTR) > 80,]
DATA$V_SEQ <- substr(DATA$SEQUENCE_INPUT, 1, (DATA$V_SEQ_START+DATA$V_SEQ_LENGTH-1))

trbv6_23_5utr_fasta <- paste0(figure_folder, "TRBV6_23.fasta")
write.fasta(as.list(DATA$UTR), DATA$SEQ_NAME, trbv6_23_5utr_fasta)


upstream_v_seqs_df <- upstream_v_seqs_df[(upstream_v_seqs_df$CLUSTER_COUNT > 2 & upstream_v_seqs_df$FRAC > 0.2) | upstream_v_seqs_df$FRAC > 0.4,]

alleles <- unique(upstream_v_seqs_df$ALLELE)

for (allele in alleles) {
  temp <- upstream_v_seqs_df[upstream_v_seqs_df$ALLELE==allele,]
  n_cons <- nrow(temp)
  if (nrow(temp) > 1) {
    for (i in 1:4){
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
UTR.db.flt <- UTR.db.flt[UTR.db.flt$GENE == "TRBV6-23",]
UTR.db.flt <- UTR.db.flt[!grepl("G47A", UTR.db.flt$ALLELE),]
UTR.db.flt$SEQ <- stri_reverse(UTR.db.flt$CONS)

write.fasta(as.list(UTR.db.flt$SEQ), UTR.db.flt$ALLELE, paste0(figure_folder, "TRBV6_23_5UTR_ref.fasta"))

# create TRBV6-23 5 UTR db reference
ref_folder <- paste0(figure_folder, "ref/")
dir.create(ref_folder)
reference_db <- paste0(ref_folder, "TRBV6_23_5UTR_ref")
system(paste0("makeblastdb -parse_seqids -dbtype nucl -in ", figure_folder, "TRBV6_23_5UTR_ref.fasta -out ", reference_db))
# run blastn
blast_output <- paste0(figure_folder, "result.tab")
system(paste0("blastn -query ", trbv6_23_5utr_fasta, " -db ", reference_db, ' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
       " -out ", blast_output))

headers <- unlist(strsplit("seq_id allele identity length mismatch gapopen qstart qend sstart send evalue bitscore", " "))

##################################################################################################################################################
###########################  Figure  #############################################################################################################
##################################################################################################################################################

TRBV6_23 <- read.delim(blast_output, header=F, sep = "\t", stringsAsFactors = F)
names(TRBV6_23) <- headers
TRBV6_23 <- TRBV6_23[TRBV6_23$identity == 100 & TRBV6_23$length > 80,]
TRBV6_23$subject <- sapply(strsplit(TRBV6_23$seq_id, "_", fixed = T), "[", 1)
TRBV6_23$g47 <- grepl("G47A", TRBV6_23$seq_id)
TRBV6_23$allele[TRBV6_23$g47] <- "TRBV6-23*01_G47A"

TRBV6_23$subject[grepl("^C[0-9]", TRBV6_23$subject)] <- paste0("H", TRBV6_23$subject[grepl("^C[0-9]", TRBV6_23$subject)])
TRBV6_23$subject <- gsub("T", "", TRBV6_23$subject)

temp <- TRBV6_23 %>% group_by(subject, allele) %>% dplyr::summarise(N=n())
total_per_s <- temp %>% group_by(subject) %>% dplyr::summarise(TOTAL=sum(N))
temp <- merge(temp, total_per_s, by = "subject")

temp$freq <- temp$N / temp$TOTAL

v623_graph <- ggplot(temp, aes(x=subject, y=freq, fill=gsub("TRBV6-23*", "", allele, fixed = T))) + 
  geom_bar(stat="identity", colour="white") + 
  scale_y_continuous(expand = c(0, 0))+
  ylab("Frequency") + labs(fill="TRBV6-2/TRBV6-3 allele") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16), axis.text = element_text(size = 14))


# TRBV7-2 figure
HCV_GENO_TRBV72 <- HCV_ALL_GENO[HCV_ALL_GENO$GENE=="TRBV7-2",]
HCV_GENO_TRBV72$SUBJECT[grepl("^C[0-9]", HCV_GENO_TRBV72$SUBJECT)] <- paste0("H", HCV_GENO_TRBV72$SUBJECT[grepl("^C[0-9]", HCV_GENO_TRBV72$SUBJECT)])
HCV_GENO_TRBV72$SUBJECT <- gsub("T", "", HCV_GENO_TRBV72$SUBJECT)

HCV_GENO_TRBV72 <- HCV_GENO_TRBV72 %>% select(SUBJECT, ALLELES, COUNTS, TOTAL) %>% separate_rows(ALLELES, COUNTS, sep=",")
HCV_GENO_TRBV72$FREQ <- as.numeric(HCV_GENO_TRBV72$COUNTS) / HCV_GENO_TRBV72$TOTAL

HCV_GENO_TRBV72 <- HCV_GENO_TRBV72[HCV_GENO_TRBV72$ALLELES == "01" | HCV_GENO_TRBV72$ALLELES == "02",]
# HCV_GENO_TRBV72$ALLELES <- paste0("V7-2*",HCV_GENO_TRBV72$ALLELES)

V72_graph <- ggplot(HCV_GENO_TRBV72, aes(x=SUBJECT, y=FREQ, fill=ALLELES)) + 
  geom_bar(stat="identity", colour="white") + 
  scale_fill_hue(l=40, c=35) +
  scale_y_continuous(expand = c(0, 0))+
  ylab("Frequency") + labs(fill="TRBV7-2 allele") + xlab("Individual") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16), axis.text = element_text(size = 14),
        axis.text.x= element_text(angle = 90, vjust = 0.5))



##################################################################################################################################################
################### Collapse the 2 plots into one figure #########################################################################################
##################################################################################################################################################

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


v623_leg <-get_legend(v623_graph)
v72_leg <-get_legend(V72_graph)

V72_graph <- V72_graph + theme(legend.position="none")
v623_graph <- v623_graph + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank())


g <- grid.arrange(v623_graph, v623_leg, V72_graph, v72_leg,
             widths = c(0.8, 0.2),
             heights = c(1, 1.2),
             layout_matrix = rbind(c(1, 2),
                                   c(3, 4)))

ggsave(paste0(figure_folder, "V6_23_V7_2_correlation.pdf"), width = 14, height = 8, plot = g)
# ggsave(paste0(figure_folder, "V6_23_V7_2_correlation.png"), width = 14, height = 8, plot = g)
