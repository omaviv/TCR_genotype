library(stringr)
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(seqinr)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_2/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

# load DS1 & DS2 genotypes mainly for panel "a."
hcv_genotypes <- read.delim(paste0(required_files_folder, "HCV_Genotypes_const_2.tab"), header = T, sep = "\t", stringsAsFactors = F)
sc_genotypes <- read.delim(paste0(required_files_folder, "Single_cell_TR_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)
biomed_genotypes <- read.delim(paste0(required_files_folder, "BIOMED2_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)
adaptive_genotypes <- read.delim(paste0(required_files_folder, "Adaptive_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)

full_novel_appearance_df <- read.delim(paste0(required_files_folder, "DS1_DS2_all_novel_occurrences.tab"), header = T, sep = "\t", stringsAsFactors = F)

TRBV_GERM <- read.fasta(paste0(references_folder, "TRBV.fasta"), as.string = TRUE)
full_trbv <- TRBV_GERM
biomed_trbv <- read.fasta(paste0(references_folder, "BIOMED2_TRBV.fasta"),as.string = TRUE)
adaptive_trbv <- read.fasta(paste0(references_folder, "Adaptive_TRBV.fasta"),as.string = TRUE)

novel_compare_df <- read.delim(paste0(required_files_folder, "full_length_undocumented_alleles.tsv"), header = T, sep = "\t", stringsAsFactors = F)
biomed_novel_df <- read.delim(paste0(required_files_folder, "bp_undocumented_alleles.tsv"), header = T, sep = "\t", stringsAsFactors = F)
adaptive_novel_df <- read.delim(paste0(required_files_folder, "ap_undocumented_alleles.tsv"), header = T, sep = "\t", stringsAsFactors = F)

# Change the names of the collapsed genes
# TRBV6-23 > TRBV6-2/TRBV6-3
GENE.loc[["TRB"]] <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", GENE.loc[["TRB"]])
hcv_genotypes$GENE <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", hcv_genotypes$GENE)
sc_genotypes$GENE <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", sc_genotypes$GENE)
biomed_genotypes$GENE <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", biomed_genotypes$GENE)
adaptive_genotypes$GENE <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", adaptive_genotypes$GENE)
full_novel_appearance_df$novel_allele <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", full_novel_appearance_df$novel_allele)
names(TRBV_GERM) <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", names(TRBV_GERM))
names(full_trbv) <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", names(full_trbv))
names(biomed_trbv) <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", names(biomed_trbv))
names(adaptive_trbv) <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", names(adaptive_trbv))
novel_compare_df$FULL_ALLELE_NAME <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", novel_compare_df$FULL_ALLELE_NAME)
biomed_novel_df$FULL_ALLELE_NAME <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", biomed_novel_df$FULL_ALLELE_NAME)
adaptive_novel_df$FULL_ALLELE_NAME <- gsub("TRBV6-23", "TRBV6-2/TRBV6-3", adaptive_novel_df$FULL_ALLELE_NAME)

# TRBV12-34 > TRBV12-3/TRBV12-4
GENE.loc[["TRB"]] <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", GENE.loc[["TRB"]])
biomed_genotypes$GENE <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", biomed_genotypes$GENE)
adaptive_genotypes$GENE <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", adaptive_genotypes$GENE)
names(biomed_trbv) <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", names(biomed_trbv))
names(adaptive_trbv) <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", names(adaptive_trbv))
biomed_novel_df$FULL_ALLELE_NAME <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", biomed_novel_df$FULL_ALLELE_NAME)
adaptive_novel_df$FULL_ALLELE_NAME <- gsub("TRBV12-34", "TRBV12-3/TRBV12-4", adaptive_novel_df$FULL_ALLELE_NAME)

# TRBV3-12 > TRBV3-1/TRBV3-2
GENE.loc[["TRB"]] <- gsub("TRBV3-12", "TRBV3-1/TRBV3-2", GENE.loc[["TRB"]])
adaptive_genotypes$GENE <- gsub("TRBV3-12", "TRBV3-1/TRBV3-2", adaptive_genotypes$GENE)
names(adaptive_trbv) <- gsub("TRBV3-12", "TRBV3-1/TRBV3-2", names(adaptive_trbv))
adaptive_novel_df$FULL_ALLELE_NAME <- gsub("TRBV3-12", "TRBV3-1/TRBV3-2", adaptive_novel_df$FULL_ALLELE_NAME)

# TRBV6-56 > TRBV6-5/TRBV6-6
GENE.loc[["TRB"]] <- gsub("TRBV6-56", "TRBV6-5/TRBV6-6", GENE.loc[["TRB"]])
adaptive_genotypes$GENE <- gsub("TRBV6-56", "TRBV6-5/TRBV6-6", adaptive_genotypes$GENE)
names(adaptive_trbv) <- gsub("TRBV6-56", "TRBV6-5/TRBV6-6", names(adaptive_trbv))
adaptive_novel_df$FULL_ALLELE_NAME <- gsub("TRBV6-56", "TRBV6-5/TRBV6-6", adaptive_novel_df$FULL_ALLELE_NAME)

########################################################################################################################


# remove repeated reads, potential sequencing error
Repeated_Read <- function(x,novel_seq){
  NT <- as.numeric(gsub("[A-Z]",'',x))
  SNPS <- gsub("[0-9]+([[:alpha:]]*).*",'\\1',x)
  SNP <- strsplit(SNPS,"")[[1]][2]
  OR_SNP <- strsplit(SNPS,"")[[1]][1]
  substr(novel_seq, NT, NT) <- OR_SNP 
  novel_seq <- c(substr(novel_seq,(NT),(NT+3)),substr(novel_seq,(NT-1),(NT+2)),
                 substr(novel_seq,(NT-2),(NT+1)),substr(novel_seq,(NT-3),(NT)))
  PAT <- paste0(c(paste0(c(rep(SNP,3),OR_SNP),collapse = ""),
                  paste0(c(rep(SNP,2),OR_SNP,SNP),collapse = ""),
                  paste0(c(SNP,OR_SNP,rep(SNP,2)),collapse = ""),
                  paste0(c(OR_SNP,rep(SNP,3)),collapse = "")), collapse = "|")
  if( any(grepl(PAT,novel_seq))) return(paste0(gsub(SNP,"X",gsub(OR_SNP,"z",novel_seq[grepl(PAT,novel_seq)])),collapse = "^"))
  else return(NA)
}

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

hcv_genotypes$SUBJECT[grepl("^C[0-9]", hcv_genotypes$SUBJECT)] <- paste0("H", hcv_genotypes$SUBJECT[grepl("^C[0-9]", hcv_genotypes$SUBJECT)])
hcv_genotypes$SUBJECT <- gsub("T", "", hcv_genotypes$SUBJECT)
hcv_genotypes$dataset <- "DS1"

sc_genotypes$dataset <- "DS2"
sc_genotypes <- sc_genotypes[substr(sc_genotypes$GENE,1,3)=="TRB",]
genotypes <- rbind(hcv_genotypes, sc_genotypes)

genotypes$SUBJECT <- paste0(genotypes$dataset, ": ", genotypes$SUBJECT)

genotypes <- genotypes[!genotypes$GENE %in% PSEUDO[["TRB"]],]

novel_appearance_df <- full_novel_appearance_df
novel_appearance_df$gene <- sapply(strsplit(novel_appearance_df$novel_allele, "*", fixed = T), "[", 1)
novel_appearance_df <- novel_appearance_df[!novel_appearance_df$novel_allele %in% PSEUDO[["TRB"]],]


k_diff_threshold = 3
y_label = "Undocumented alleles"

# read the merged genotype file
NOVEL_ALLELE_DF <- genotypes
# Filter as you wish
# NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[NOVEL_ALLELE_DF$K_DIFF >= k_diff_threshold,]
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[grepl("_[A-Z]", NOVEL_ALLELE_DF$GENOTYPED_ALLELES),]

# Extracts the novel alleles
NOVEL_ALLELE_DF$HETERO <- grepl(",", NOVEL_ALLELE_DF$GENOTYPED_ALLELES)
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF %>% separate_rows(GENOTYPED_ALLELES, sep = ",")
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[grepl("_[A-Z]", NOVEL_ALLELE_DF$GENOTYPED_ALLELES),]
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF %>% separate_rows(ALLELES, COUNTS, sep = ",")
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[NOVEL_ALLELE_DF$ALLELES == NOVEL_ALLELE_DF$GENOTYPED_ALLELES,]
NOVEL_ALLELE_DF$FULL_ALLELE_NAME <- paste0(NOVEL_ALLELE_DF$GENE, "*", NOVEL_ALLELE_DF$GENOTYPED_ALLELES)


# allele notes
# names(TRBV_GERM) <- sapply(strsplit(names(TRBV_GERM),"|",fixed = T),"[",2) 
TRBV_GERM <- toupper(unlist(TRBV_GERM))

unique_novel_alleles <- unique(NOVEL_ALLELE_DF$FULL_ALLELE_NAME)
for (novel in unique_novel_alleles) {
  known_allele <- sapply(strsplit(novel, "_", fixed = T), "[", 1)
  novel_seq <- TRBV_GERM[[known_allele]]
  snps <- strsplit(novel, "_", fixed = T)[[1]]
  snps <- snps[-1]
  for (snp in snps) {
    pos <- as.numeric(str_extract(snp, "[0-9]+"))
    new_nt <- str_extract(snp, "[A-Z]$")
    substr(novel_seq, pos, pos) <- new_nt
  }
  TRBV_GERM[novel] <- novel_seq
}


# write.fasta(sequences = as.list(c(TRBV_GERM)), names = names(TRBV_GERM), "TRBV_with_novel.fasta", open="w")
# 
# write.fasta(sequences=as.list(gsub(TRBV_GERM, pattern = '.', replacement = '', fixed = T)),
#             names=names(TRBV_GERM), "TRBV_with_novel_ref.fasta", open="w")


novel_checks <- data.frame(ALLELE = unique_novel_alleles, SNP_XXXX = "", stringsAsFactors = F)

print("check RR")
novel_checks$SNP_XXXX <- unlist(sapply(1:nrow(novel_checks), function(i){
  subs <- strsplit(novel_checks$ALLELE[i],'_')[[1]][-1]
  allele <- strsplit(novel_checks$ALLELE[i],'_')[[1]][1]
  novel_seq <- TRBV_GERM[[allele]]
  for (snp in subs) {
    pos <- as.numeric(str_extract(snp, "[0-9]+"))
    new_nt <- str_extract(snp, "[A-Z]$")
    substr(novel_seq, pos, pos) <- new_nt
  }
  RR <- sapply(subs,Repeated_Read,novel_seq=novel_seq,simplify = F)
  if(sum(!is.na(RR))!=0) paste0("Repeated Read: ", paste0(names(RR),"-",paste0(RR,collapse = '^'),collapse = "/"),";")
  else return("")
}))
novel_checks$padded <- ifelse(grepl("R", novel_checks$SNP_XXXX), T, F)

print("Check unary model")
novel_checks$unary <- unlist(lapply(unique_novel_alleles, function(allele) { 
  (max(novel_appearance_df$allele_freq[novel_appearance_df$novel_allele==allele & !is.na(novel_appearance_df$allele_freq) & novel_appearance_df$carry != "Not carry"]) < 0.2) & 
    (nrow(novel_appearance_df[novel_appearance_df$allele_freq >0,]) > 2)
}))
novel_checks$both <- novel_checks$padded & novel_checks$unary


# Extracts the subjects and the novel alleles
subjects <- unique(NOVEL_ALLELE_DF$SUBJECT)
novel_alleles <- unique(NOVEL_ALLELE_DF$FULL_ALLELE_NAME)

# Adapt the table for the heatmap graph
SUBJECT_NOVEL_ALLELES_DF <- data.frame()
subject_number <- length(subjects)
for (novel_allele in novel_alleles) {
  temp_df <- data.frame("SUBJECT" = subjects,"NOVEL_ALLELES" = rep(novel_allele, subject_number), "HAVING" = rep(0, subject_number))
  temp_df$HAVING[temp_df$SUBJECT %in% NOVEL_ALLELE_DF$SUBJECT[NOVEL_ALLELE_DF$FULL_ALLELE_NAME==novel_allele & NOVEL_ALLELE_DF$HETERO]] <- 1
  temp_df$HAVING[temp_df$SUBJECT %in% NOVEL_ALLELE_DF$SUBJECT[NOVEL_ALLELE_DF$FULL_ALLELE_NAME==novel_allele & (!NOVEL_ALLELE_DF$HETERO)]] <- 2
  
  SUBJECT_NOVEL_ALLELES_DF <- rbind(SUBJECT_NOVEL_ALLELES_DF, temp_df)
}

# Order the subjects and the novel alleles
NOVEL_DIS_DF <- NOVEL_ALLELE_DF %>% group_by(GENE, FULL_ALLELE_NAME) %>% dplyr::summarise(SEQUENCE_NUM = sum(as.integer(COUNTS)), SUBJECTS_NUM = n())
NOVEL_DIS_DF$HETEROZYGOUS <- unlist(lapply(NOVEL_DIS_DF$FULL_ALLELE_NAME, function(x) {
  nrow(NOVEL_ALLELE_DF[NOVEL_ALLELE_DF$FULL_ALLELE_NAME == x & NOVEL_ALLELE_DF$HETERO,])
}))

NOVEL_DIS_DF$HOMOZYGOUS <- NOVEL_DIS_DF$SUBJECTS_NUM - NOVEL_DIS_DF$HETEROZYGOUS

allele_order <- NOVEL_DIS_DF$FULL_ALLELE_NAME[order(-NOVEL_DIS_DF$SUBJECTS_NUM)]

subject_order <- lapply(subjects, function(subject) {
  sum(unlist(lapply(1:length(allele_order), function(i) {
    SUBJECT_NOVEL_ALLELES_DF$HAVING[SUBJECT_NOVEL_ALLELES_DF$SUBJECT == subject & SUBJECT_NOVEL_ALLELES_DF$NOVEL_ALLELES == allele_order[i]] * (3 ^ (length(allele_order) - i))
  })))
})
names(subject_order) <- subjects
subject_order <- subject_order[order(-unlist(subject_order))]
subject_order <- names(subject_order)

# converting the list values to charecter
SUBJECT_NOVEL_ALLELES_DF$HAVING[SUBJECT_NOVEL_ALLELES_DF$HAVING == 0] <- ""
SUBJECT_NOVEL_ALLELES_DF$HAVING[SUBJECT_NOVEL_ALLELES_DF$HAVING == 1] <- "Heterozygous"
SUBJECT_NOVEL_ALLELES_DF$HAVING[SUBJECT_NOVEL_ALLELES_DF$HAVING == 2] <- "Homozygous"

SUBJECT_NOVEL_ALLELES_DF$HAVING <- as.character(SUBJECT_NOVEL_ALLELES_DF$HAVING)
SUBJECT_NOVEL_ALLELES_DF <- SUBJECT_NOVEL_ALLELES_DF[order(factor(SUBJECT_NOVEL_ALLELES_DF$NOVEL_ALLELES, levels = allele_order),
                                                           factor(SUBJECT_NOVEL_ALLELES_DF$SUBJECT, levels = subject_order)),]

dataset_order <- unlist(lapply(subject_order, function(subject){
  genotypes$dataset[genotypes$SUBJECT == subject][[1]]
}))
subjects_colors <- ifelse(dataset_order=="DS1", "blue", "red")
print(subjects_colors)

allele_colors <- unlist(lapply(allele_order, function(allele) {
  if (novel_checks$both[novel_checks$ALLELE==allele]) {
    return("navy")
  } else if (novel_checks$padded[novel_checks$ALLELE==allele]) {
    return("green4")
  } else if (novel_checks$unary[novel_checks$ALLELE==allele]) {
    return("red")
  }
  return("black")
}))


# Generate the graph
novel_allele_graph <- ggplot(SUBJECT_NOVEL_ALLELES_DF, aes(x = SUBJECT, y = NOVEL_ALLELES, fill = HAVING, group = HAVING)) +
  geom_tile() + 
  scale_x_discrete(limits = subject_order) +
  scale_y_discrete(limits = allele_order) +
  scale_fill_manual(values = c("white", "skyblue", "blue")) +
  xlab("Individuals") + ylab(y_label) +
  guides(fill=guide_legend(title="", reverse=TRUE)) + 
  theme(axis.title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=18),
        axis.text.x= element_text(size=14, angle = 90, vjust = 0.3),
        axis.text.y= element_text(size=14, color = allele_colors))
# axis.text.x= element_text(angle = 90, vjust = 0.3))
novel_allele_graph


lower_than_k <- NOVEL_ALLELE_DF[NOVEL_ALLELE_DF$K_DIFF<k_diff_threshold,]
lines_df <- c()
for (i in 1:nrow(lower_than_k)) {
  subject <- lower_than_k$SUBJECT[i]
  allele <- lower_than_k$FULL_ALLELE_NAME[i]
  x_ind  <- which(subject_order==subject)
  y_ind <- which(allele_order==allele)
  n_sides <- 3
  
  gaps <- c(1:(n_sides)/(n_sides+1), rep(0, n_sides+1))
  x_starts <- ((x_ind-0.5)+gaps)
  y_starts <- ((y_ind-0.5)+rev(gaps))
  y_ends <- ((y_ind+0.5)-gaps)
  x_ends <- ((x_ind+0.5)-rev(gaps))
  square_lines_df <- data.frame(x=x_starts, y=y_starts, xend=x_ends, yend=y_ends)
  lines_df <- rbind(lines_df,square_lines_df)
}

novel_allele_graph <- novel_allele_graph + annotate("segment", 
                                                    x = lines_df$x,
                                                    xend = lines_df$xend,
                                                    y = lines_df$y,
                                                    yend = lines_df$yend,
                                                    color = "gray50",
                                                    size=0.3)

plot(novel_allele_graph)

##########################################################################################################################

# novel_appearance_df <- novel_appearance_df[novel_appearance_df$gene_single_allele_assignments > 30,]
novel_appearance_df <- full_novel_appearance_df
novel_appearance_df$n_sequences <- unlist(lapply(novel_appearance_df$gene_single_allele_assignments, function(x){
  if (x < 10) {
    return("<10")
  } else if (x <= 30) {
    return("10-30")
  } else {
    return("30<")
  }
}))
novel_appearance_df$n_sequences <- factor(novel_appearance_df$n_sequences, levels = c("<10", "10-30", "30<"))

snp_freq_plot <- ggplot(novel_appearance_df, aes(y=novel_allele, x=allele_freq, color = carry, shape=n_sequences)) + 
  geom_point() + theme_bw()+
  scale_y_discrete(limits = allele_order) +
  scale_color_manual(values = c("skyblue", "blue", "grey")) +
  xlab('Allele fraction') + ylab('Polymorphism call') + 
  geom_hline(yintercept = 0.25, linetype="dashed", color = "gray50") +
  theme(axis.text.x=element_text(size=14),
    axis.title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=18))
snp_freq_plot


shape_legend <- ggplot(novel_appearance_df, aes(y=novel_allele, x=allele_freq, shape=n_sequences)) + 
  geom_point(size=3) + labs(shape="Number of gene single\nassignment sequences") + 
  theme(axis.title = element_text(size=20), legend.title = element_text(size=18), legend.text = element_text(size=16))

shape_legend <- get_legend(shape_legend)

snp_freq_plot <- snp_freq_plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
temp <- novel_allele_graph  + theme(legend.position = "none")

# build legend
leg_df <- data.frame(x=rep(c(1:4), 4), y=unlist(lapply(1:4, rep, 4)), stringsAsFactors = F)
leg_df$carry <- ""
leg_texts <- c("", "lk < 3", "Undocumented allele\nnot in genotype", "Heterozygous to\nthe undocumented allele", "Homozygous to\nthe undocumented allele")
leg_df$carry[leg_df$x == 1] <- leg_texts[-1]
leg_df$carry <- factor(leg_df$carry, levels = leg_texts)

leg_plot <- ggplot(leg_df, aes(x = x, y = y, fill = carry, label=carry)) +
  geom_tile() + 
  geom_text(nudge_x = 3.2, size=5.4) + 
  scale_fill_manual(values = c("white", "white", "grey", "skyblue", "blue")) +
  theme_nothing()

n_sides <- 3
gaps <- c(1:(n_sides)/(n_sides+1), rep(0, n_sides+1))
x_starts <- ((0.5)+gaps)
y_starts <- ((0.5)+rev(gaps))
y_ends <- ((1.5)-gaps)
x_ends <- ((1.5)-rev(gaps))

leg_plot <- leg_plot + annotate("segment", 
                                x = x_starts,
                                xend = x_ends,
                                y = y_starts,
                                yend = y_ends,
                                color = "gray50",
                                size=0.3)


# grid_legends <- grid.arrange(geno_leg, snp_freq_plot_led, nrow=2)

p1.common.y <- ggplot_gtable(ggplot_build(temp))
p2.common.y <- ggplot_gtable(ggplot_build(snp_freq_plot))

# copy the plot height from p1 to p2
p2.common.y$heights <- p1.common.y$heights

ds1_ds2_graph <- grid.arrange(p1.common.y, p2.common.y, ncol=2, widths=c(0.75,0.25))

ds1_ds2_graph <- grid.arrange(grobs = list(ds1_ds2_graph, leg_plot, shape_legend), layout_matrix = rbind(c(1,1,1,1,1,3),
                                                                                 c(1,1,1,1,1,2),
                                                                                 c(1,1,1,1,1,NA)))



#####################################################################################################################
############################ Allele appearance in DS1&DS2 #################################################################
#####################################################################################################################
full_trbv <- as.character(names(full_trbv))
full_trbv <- sapply(base::strsplit(full_trbv, "*", fixed = T), "[", 1)
full_trbv <- data.frame(table(full_trbv), stringsAsFactors = F)
names(full_trbv) <- c("GENE", "allele_number")
full_trbv$GENE <- as.character(full_trbv$GENE)
full_trbv$novel <- "Documented alleles"
full_trbv$group <- "IMGT"

gene2show <- GENE.loc[["TRB"]][!GENE.loc[["TRB"]] %in% unique(full_trbv$GENE)]
gene2show <- gene2show[grepl("V", gene2show)]
gene2show <- gene2show[!gene2show %in% PSEUDO[["TRB"]]]
for (gene in gene2show) {
  full_trbv[nrow(full_trbv)+1,] <- c(gene, 0, "Documented alleles", "IMGT")
}

hcv_genotypes <- hcv_genotypes[substr(hcv_genotypes$GENE,1,4)=="TRBV",]
sc_genotypes <- sc_genotypes[substr(sc_genotypes$GENE,1,4)=="TRBV",]

genotypes <- rbind(hcv_genotypes, sc_genotypes)
genotypes <- genotypes[!grepl("Del", genotypes$GENOTYPED_ALLELES),]
genotypes <- genotypes %>% separate_rows(GENOTYPED_ALLELES, sep = ",")
genotypes$GENOTYPED_ALLELES[grepl("_[0-9]",genotypes$GENOTYPED_ALLELES)] <- 
  sapply(strsplit(genotypes$GENOTYPED_ALLELES[grepl("_[0-9]",genotypes$GENOTYPED_ALLELES)], "_", fixed=T), "[", 1)

novel2drop <- novel_compare_df$FULL_ALLELE_NAME[novel_compare_df$padded | novel_compare_df$unary]
for (n2d in novel2drop) {
  gene <- sapply(strsplit(n2d, "*", fixed = T), "[", 1) 
  allele <- sapply(strsplit(n2d, "*", fixed = T), "[", 2) 
  genotypes <- genotypes[!(genotypes$GENE==gene & genotypes$GENOTYPED_ALLELES==allele),]
}

allele_appearence <- genotypes %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N=n())
allele_appearence$novel <- grepl("_", allele_appearence$GENOTYPED_ALLELES)
allele_appearence <- allele_appearence %>% group_by(GENE, novel) %>% dplyr::summarise(allele_number=n())
allele_appearence$allele_number[allele_appearence$novel] <- unlist(lapply(allele_appearence$GENE[allele_appearence$novel], 
                                                                          function(gene){
                                                                            sum(allele_appearence$allele_number[allele_appearence$GENE==gene])
                                                                          }))
allele_appearence$novel <- ifelse(allele_appearence$novel, "Undocumented alleles", "Observed documented alleles")
allele_appearence$group <- "genotype"


full_trbv <- rbind(full_trbv, as.data.frame(allele_appearence))

full_trbv$novel <- factor(full_trbv$novel, c("Documented alleles", "Observed documented alleles", "Undocumented alleles"))
full_trbv$group <- factor(full_trbv$group, c("IMGT", "genotype"))

full_trbv <- full_trbv[!full_trbv$GENE %in% PSEUDO[["TRB"]],] # throws pseudo genes
full_trbv <- full_trbv[!grepl("OR", full_trbv$GENE),] # throws orphan genes
full_trbv$GENE <- gsub("TRB", "", full_trbv$GENE)
full_trbv$GENE <- factor(full_trbv$GENE, levels = gsub("TRB", "", GENE.loc[["TRB"]]))

full_trbv$allele_number <- as.integer(full_trbv$allele_number)
full_trbv <- full_trbv[!full_trbv$GENE %in% PSEUDO[["TRB"]],]
full_gene_allele_num_plot <- ggplot(data=full_trbv[order(full_trbv$allele_number, decreasing = T),], aes(x=GENE, y=allele_number, fill=novel, group=group)) +
  geom_bar(width = 0.8, stat="identity", color="black", position=position_dodge())+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(full_trbv$allele_number))) +
  scale_fill_brewer(palette = "Dark2") + 
  xlab("TRBV Genes") + ylab("Number of alleles") + ggtitle("DS1 & DS2") + 
  guides(fill=guide_legend("", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=14, angle = 90, vjust = 0.5), axis.text.y = element_text(size=14),
        axis.title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=18),
        legend.position=c(0.5,0.9))

full_gene_allele_num_plot

#####################################################################################################################
############################ Allele appearance in DS3 #################################################################
#####################################################################################################################

biomed_trbv <- as.character(names(biomed_trbv))
biomed_trbv <- sapply(base::strsplit(biomed_trbv, "*", fixed = T), "[", 1)
biomed_trbv <- as.data.frame(table(biomed_trbv))
names(biomed_trbv) <- c("GENE", "allele_number")
biomed_trbv$GENE <- as.character(biomed_trbv$GENE)
biomed_trbv$novel <- "Documented alleles"
biomed_trbv$group <- "IMGT"

gene2show <- GENE.loc[["TRB"]][!GENE.loc[["TRB"]] %in% unique(biomed_trbv$GENE)]
gene2show <- gene2show[grepl("V", gene2show)]
gene2show <- gene2show[!gene2show %in% PSEUDO[["TRB"]]]
for (gene in gene2show) {
  biomed_trbv[nrow(biomed_trbv)+1,] <- c(gene, 0, "Documented alleles", "IMGT")
}


biomed_genotypes <- biomed_genotypes[substr(biomed_genotypes$GENE,1,4)=="TRBV",]
biomed_genotypes <- biomed_genotypes[!grepl("Del", biomed_genotypes$GENOTYPED_ALLELES),]
biomed_genotypes <- biomed_genotypes %>% separate_rows(GENOTYPED_ALLELES, sep = ",")
biomed_genotypes$GENOTYPED_ALLELES[grepl("_[0-9]",biomed_genotypes$GENOTYPED_ALLELES)] <- 
  sapply(strsplit(biomed_genotypes$GENOTYPED_ALLELES[grepl("_[0-9]",biomed_genotypes$GENOTYPED_ALLELES)], "_", fixed=T), "[", 1)

novel2drop <- biomed_novel_df$FULL_ALLELE_NAME[biomed_novel_df$padded | biomed_novel_df$unary]
for (n2d in novel2drop) {
  gene <- sapply(strsplit(n2d, "*", fixed = T), "[", 1) 
  allele <- sapply(strsplit(n2d, "*", fixed = T), "[", 2) 
  biomed_genotypes <- biomed_genotypes[!(biomed_genotypes$GENE==gene & biomed_genotypes$GENOTYPED_ALLELES==allele),]
}


allele_appearence <- biomed_genotypes %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N=n())
allele_appearence$novel <- grepl("_", allele_appearence$GENOTYPED_ALLELES)
allele_appearence <- allele_appearence %>% group_by(GENE, novel) %>% dplyr::summarise(allele_number=n())
allele_appearence$allele_number[allele_appearence$novel] <- unlist(lapply(allele_appearence$GENE[allele_appearence$novel], 
                                                                          function(gene){
                                                                            sum(allele_appearence$allele_number[allele_appearence$GENE==gene])
                                                                          }))
allele_appearence$novel <- ifelse(allele_appearence$novel, "Undocumented alleles", "Observed documented alleles")
allele_appearence$group <- "genotype"

biomed_trbv <- rbind(biomed_trbv, as.data.frame(allele_appearence))

biomed_trbv$novel <- factor(biomed_trbv$novel, c("Documented alleles", "Observed documented alleles", "Undocumented alleles"))
biomed_trbv$group <- factor(biomed_trbv$group, c("IMGT", "genotype"))

biomed_trbv <- biomed_trbv[!biomed_trbv$GENE %in% PSEUDO[["TRB"]],] # throws pseudo genes
biomed_trbv <- biomed_trbv[!grepl("OR", biomed_trbv$GENE),] # throws orphan genes
biomed_trbv$GENE <- gsub("TRB", "", biomed_trbv$GENE)
biomed_trbv$GENE <- factor(biomed_trbv$GENE, levels = gsub("TRB", "", GENE.loc[["TRB"]]))

biomed_trbv <- biomed_trbv[!biomed_trbv$GENE %in% PSEUDO[["TRB"]],]
biomed_trbv$allele_number <- as.integer(biomed_trbv$allele_number)
biomed_gene_allele_num_plot <- ggplot(data=biomed_trbv[order(biomed_trbv$allele_number, decreasing = T),], aes(x=GENE, y=allele_number, fill=novel, group=group)) +
  geom_bar(width = 0.8, stat="identity", color="black", position=position_dodge())+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Dark2") + 
  xlab("TRBV Genes") + ylab("Number of alleles") + ggtitle("DS3") + 
  guides(fill=guide_legend("", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=14, angle = 90, vjust = 0.5), axis.text.y = element_text(size=14),
        axis.title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=18),
        legend.position=c(0.5,0.9))

biomed_gene_allele_num_plot

#####################################################################################################################
############################ Allele appearance in DS4 #################################################################
#####################################################################################################################

adaptive_trbv <- as.character(names(adaptive_trbv))
adaptive_trbv <- sapply(base::strsplit(adaptive_trbv, "*", fixed = T), "[", 1)
adaptive_trbv <- as.data.frame(table(adaptive_trbv))
names(adaptive_trbv) <- c("GENE", "allele_number")
adaptive_trbv$GENE <- as.character(adaptive_trbv$GENE)
adaptive_trbv$novel <- "Documented alleles"
adaptive_trbv$group <- "IMGT"

gene2show <- GENE.loc[["TRB"]][!GENE.loc[["TRB"]] %in% unique(adaptive_trbv$GENE)]
gene2show <- gene2show[grepl("V", gene2show)]
gene2show <- gene2show[!gene2show %in% PSEUDO[["TRB"]]]
for (gene in gene2show) {
  adaptive_trbv[nrow(adaptive_trbv)+1,] <- c(gene, 0, "Documented alleles", "IMGT")
}

adaptive_genotypes <- adaptive_genotypes[substr(adaptive_genotypes$GENE,1,4)=="TRBV",]
adaptive_genotypes <- adaptive_genotypes[!grepl("Del", adaptive_genotypes$GENOTYPED_ALLELES),]
adaptive_genotypes <- adaptive_genotypes %>% separate_rows(GENOTYPED_ALLELES, sep = ",")

novel2drop <- adaptive_novel_df$FULL_ALLELE_NAME[adaptive_novel_df$padded | adaptive_novel_df$unary]
for (n2d in novel2drop) {
  gene <- sapply(strsplit(n2d, "*", fixed = T), "[", 1) 
  allele <- sapply(strsplit(n2d, "*", fixed = T), "[", 2) 
  adaptive_genotypes <- adaptive_genotypes[!(adaptive_genotypes$GENE==gene & adaptive_genotypes$GENOTYPED_ALLELES==allele),]
}


allele_appearence <- adaptive_genotypes %>% group_by(GENE, GENOTYPED_ALLELES) %>% dplyr::summarise(N=n())
allele_appearence$novel <- grepl("_", allele_appearence$GENOTYPED_ALLELES)
allele_appearence <- allele_appearence %>% group_by(GENE, novel) %>% dplyr::summarise(allele_number=n())
allele_appearence$allele_number[allele_appearence$novel] <- unlist(lapply(allele_appearence$GENE[allele_appearence$novel], 
                                                                          function(gene){
                                                                            sum(allele_appearence$allele_number[allele_appearence$GENE==gene])
                                                                          }))
allele_appearence$novel <- ifelse(allele_appearence$novel, "Undocumented alleles", "Observed documented alleles")
allele_appearence$group <- "genotype"

adaptive_trbv <- rbind(adaptive_trbv, as.data.frame(allele_appearence))

adaptive_trbv$novel <- factor(adaptive_trbv$novel, c("Documented alleles", "Observed documented alleles", "Undocumented alleles"))
adaptive_trbv$group <- factor(adaptive_trbv$group, c("IMGT", "genotype"))

adaptive_trbv <- adaptive_trbv[!adaptive_trbv$GENE %in% PSEUDO[["TRB"]],] # throws pseudo genes
adaptive_trbv <- adaptive_trbv[!grepl("OR", adaptive_trbv$GENE),] # throws orphan genes
adaptive_trbv$GENE <- gsub("TRB", "", adaptive_trbv$GENE)
adaptive_trbv$GENE <- factor(adaptive_trbv$GENE, levels = gsub("TRB", "", GENE.loc[["TRB"]]))

adaptive_trbv$allele_number <- as.integer(adaptive_trbv$allele_number)
adaptive_gene_allele_num_plot <- ggplot(data=adaptive_trbv[order(adaptive_trbv$allele_number, decreasing = T),], aes(x=GENE, y=allele_number, fill=novel, group=group)) +
  geom_bar(width = 0.8, stat="identity", color="black", position=position_dodge())+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Dark2") + 
  xlab("TRBV Genes") + ylab("Number of alleles") + ggtitle("DS4") + 
  guides(fill=guide_legend("", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=14, angle = 90, vjust = 0.5), axis.text.y = element_text(size=14),
        axis.title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=18),
        legend.position=c(0.5,0.9))

adaptive_gene_allele_num_plot

#####################################################################################################################
############################ Group to one figure #################################################################
#####################################################################################################################

require(grid)
a_label <- textGrob(
  label = "A.",
  gp = gpar(fontsize = 34), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

b_label <- textGrob(
  label = "B.",
  gp = gpar(fontsize = 34), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

y_label <- textGrob(
  label = "Number of Alleles",
  gp = gpar(fontsize = 20), rot=90)


ds1_ds2_graph <- arrangeGrob(ds1_ds2_graph, top = a_label)

gene_allele_num_leg <- get_legend(full_gene_allele_num_plot)
full_gene_allele_num_plot <- full_gene_allele_num_plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")
biomed_gene_allele_num_plot <- biomed_gene_allele_num_plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")
adaptive_gene_allele_num_plot <- adaptive_gene_allele_num_plot + theme(axis.title.y = element_blank(), legend.position = "none")

gene_allele_num_plot <- grid.arrange(full_gene_allele_num_plot, biomed_gene_allele_num_plot, adaptive_gene_allele_num_plot, ncol=1, heights=c(1.1,1,1.8))
gene_allele_num_plot <- arrangeGrob(gene_allele_num_plot, top = b_label, left = y_label, bottom = textGrob(label = "TRBV genes", gp = gpar(fontsize = 20)))
gene_allele_num_plot <- grid.arrange(gene_allele_num_plot, gene_allele_num_leg, ncol=1, heights=c(10, 1))

g <- grid.arrange(
  ds1_ds2_graph, gene_allele_num_plot,
  heights=c(1,1),
  layout_matrix = rbind(c(1),
                        c(2))
)

ggsave(paste0(figure_folder, "Allele_appearance.pdf"), plot = g, width = 20, height = 18)
ggsave(paste0(figure_folder, "Allele_appearance.png"), plot = g, width = 15, height = 15)




