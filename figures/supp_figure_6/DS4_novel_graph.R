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
figure_folder <- paste0(project_folder, "figures/supp_figure_6/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

genotypes <- read.delim(paste0(required_files_folder, "Adaptive_All_Genotypes.tab"), header = T, sep = "\t", stringsAsFactors = F)
novel_appearance_df <- read.delim(paste0(required_files_folder, "DS4_all_novel_occurrences.tab"), header = T, sep = "\t", stringsAsFactors = F)
TRBV_imgt_ref <- paste0(references_folder, "Adaptive_TRBV.fasta")

extra_red_alleles <- c("TRBV7-4*ap01_G291C_A297G", "TRBV7-4*ap01_G291C_A297G_C314T", "TRBV7-9*ap01_G313T")
#####################################################################################################################
###########################  Figure  ################################################################################
#####################################################################################################################

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



novel_appearance_df$gene <- sapply(strsplit(novel_appearance_df$novel_allele, "*", fixed = T), "[", 1)
novel_appearance_df <- novel_appearance_df[!novel_appearance_df$gene %in% PSEUDO[["TRB"]],]

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
shape_legend <- ggplot(novel_appearance_df, aes(y=novel_allele, x=allele_freq, shape=n_sequences)) + 
  geom_point(size=3) + labs(shape="Number of gene single assignment sequences")
shape_legend <- get_legend(shape_legend)


genotypes <- genotypes[!genotypes$GENE %in% PSEUDO[["TRB"]],]

k_diff_threshold = 3

# Filter as you wish
NOVEL_ALLELE_DF <- genotypes
# NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[NOVEL_ALLELE_DF$K_DIFF >= 3,]
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[grepl("_[A-Z]", NOVEL_ALLELE_DF$GENOTYPED_ALLELES),]

# Extracts the novel alleles
NOVEL_ALLELE_DF$HETERO <- grepl(",", NOVEL_ALLELE_DF$GENOTYPED_ALLELES)
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF %>% separate_rows(GENOTYPED_ALLELES, sep = ",")
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[grepl("_", NOVEL_ALLELE_DF$GENOTYPED_ALLELES),]
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF %>% separate_rows(ALLELES, COUNTS, sep = ",")
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[NOVEL_ALLELE_DF$ALLELES == NOVEL_ALLELE_DF$GENOTYPED_ALLELES,]
NOVEL_ALLELE_DF$FULL_ALLELE_NAME <- paste0(NOVEL_ALLELE_DF$GENE, "*", NOVEL_ALLELE_DF$GENOTYPED_ALLELES)
NOVEL_ALLELE_DF$low_kdiff <- NOVEL_ALLELE_DF$K_DIFF <= k_diff_threshold

# allele notes
TRBV_GERM <- read.fasta(TRBV_imgt_ref, as.string = TRUE)
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

cutted_trbv <- substr(TRBV_GERM, 1, 316)
existing_alleles <- c()
for (allele in cutted_trbv[grepl("_", names(cutted_trbv))]) {
  if (length(cutted_trbv[cutted_trbv==allele]) > 1) {
    existing_alleles <- c(existing_alleles, cutted_trbv[cutted_trbv==allele])
    print(names(cutted_trbv[cutted_trbv==allele]))
  }
}

allele2drop <- names(existing_alleles)
allele2drop <- unique(allele2drop[grepl("_", allele2drop)])
NOVEL_ALLELE_DF <- NOVEL_ALLELE_DF[!NOVEL_ALLELE_DF$FULL_ALLELE_NAME %in% allele2drop,]
TRBV_GERM <- TRBV_GERM[!names(TRBV_GERM) %in% allele2drop]

# write.fasta(sequences = as.list(c(TRBV_GERM)), names = names(TRBV_GERM),
#             paste0(figure_folder, "Adaptive_TRBV_with_novel.fasta"), open="w")
# 
# write.fasta(sequences=as.list(gsub(TRBV_GERM, pattern = '.', replacement = '', fixed = T)),
#             names=names(TRBV_GERM),
#             paste0(figure_folder, "Adaptive_TRBV_with_novel_ref.fasta"), open="w")


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

# Order the subjects and the novel alleles
save_novel_geno <- NOVEL_ALLELE_DF
NOVEL_DIS_DF <- NOVEL_ALLELE_DF %>% group_by(GENE, FULL_ALLELE_NAME) %>% 
  dplyr::summarise(SEQUENCE_NUM = sum(as.integer(COUNTS)), SUBJECTS_NUM = n(), HETEROZYGOUS=sum(HETERO), HOMOZYGOUS=sum(!HETERO),
                   HETEROZYGOUS_LOW_K=sum(HETERO&low_kdiff), HOMOZYGOUS_LOW_K=sum((!HETERO)&low_kdiff))

high_freq_novel_df <- NOVEL_DIS_DF[NOVEL_DIS_DF$SUBJECTS_NUM >= 50,]
high_freq_allele_order <- high_freq_novel_df$FULL_ALLELE_NAME[order(-high_freq_novel_df$SUBJECTS_NUM)]

low_freq_novel_df <- NOVEL_DIS_DF[NOVEL_DIS_DF$SUBJECTS_NUM < 50,]
low_freq_allele_order <- low_freq_novel_df$FULL_ALLELE_NAME[order(-low_freq_novel_df$SUBJECTS_NUM)]


################################################## GRAPH FOR LOW ######################################################################################

# split the dataframe into 2 groups homozygous to the novel allele and heterozygous to the novel allele
homozygous_novel_df <- low_freq_novel_df
homozygous_novel_df$SUBJECTS_NUMBER <- homozygous_novel_df$HOMOZYGOUS
homozygous_novel_df$GROUP <- "Homozygous"
homozygous_novel_df <- homozygous_novel_df %>% select(GENE, FULL_ALLELE_NAME, SUBJECTS_NUMBER, GROUP)

heterozygous_novel_df <- low_freq_novel_df
heterozygous_novel_df$SUBJECTS_NUMBER <- heterozygous_novel_df$HETEROZYGOUS
heterozygous_novel_df$GROUP <- "Heterozygous"
heterozygous_novel_df <- heterozygous_novel_df %>% select(GENE, FULL_ALLELE_NAME, SUBJECTS_NUMBER, GROUP)

NOVEL_ALLELE_DF <- rbind(homozygous_novel_df, heterozygous_novel_df)
NOVEL_ALLELE_DF$FULL_ALLELE_NAME <- factor(NOVEL_ALLELE_DF$FULL_ALLELE_NAME, low_freq_allele_order)

colors <- unlist(lapply(low_freq_allele_order, function(allele) {
  if (novel_checks$both[novel_checks$ALLELE==allele]) {
    return("navy")
  } else if (novel_checks$padded[novel_checks$ALLELE==allele]) {
    return("green4")
  } else if (novel_checks$unary[novel_checks$ALLELE==allele]) {
    return("red")
  }
  return("black")
}))

low_novel_allele_graph <- ggplot(NOVEL_ALLELE_DF, aes(fill=GROUP, y=FULL_ALLELE_NAME, x=SUBJECTS_NUMBER)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("skyblue", "blue")) +
  scale_x_continuous(expand = c(0, 0), breaks = (0:4)*5) +
  guides(fill=guide_legend(title="", reverse=TRUE)) +
  ylab("Undocumented pattern alleles")  + xlab("Number of individuals") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text.y= element_text(color = colors))

low_novel_allele_graph


lines_df <- c()
NOVEL_DIS_DF <- low_freq_novel_df
for (i in 1:nrow(NOVEL_DIS_DF)) {
  allele <- NOVEL_DIS_DF$FULL_ALLELE_NAME[i]
  y_ind <- which(low_freq_allele_order==allele)
  n_sides <- 3
  
  y_gaps <- c((1:(n_sides)/(n_sides+1))*0.84, rep(0, n_sides+1))
  y_starts <- ((y_ind-0.42)+rev(y_gaps))
  y_ends <- ((y_ind+0.42)-y_gaps)
  
  n_homozygous <- NOVEL_DIS_DF$HOMOZYGOUS[i]
  n_hetero <- NOVEL_DIS_DF$HETEROZYGOUS[i]
  
  if (NOVEL_DIS_DF$HOMOZYGOUS_LOW_K[i] > 0) {
    n_low_k <- NOVEL_DIS_DF$HOMOZYGOUS_LOW_K[i]
    x_gaps <- c((1:(n_sides)/(n_sides+1))*n_low_k, rep(0, n_sides+1))
    x_starts <- ((n_homozygous-n_low_k)+x_gaps)
    x_ends <- ((n_homozygous)-rev(x_gaps))
    square_lines_df <- data.frame(x=x_starts, y=y_starts, xend=x_ends, yend=y_ends)
    lines_df <- rbind(lines_df,square_lines_df)
  }
  
  if (NOVEL_DIS_DF$HETEROZYGOUS_LOW_K[i] > 0) {
    n_low_k <- NOVEL_DIS_DF$HETEROZYGOUS_LOW_K[i]
    x_gaps <- c((1:(n_sides)/(n_sides+1))*n_low_k, rep(0, n_sides+1))
    x_starts <- ((n_homozygous+n_hetero-n_low_k)+x_gaps)
    x_ends <- ((n_homozygous+n_hetero)-rev(x_gaps))
    square_lines_df <- data.frame(x=x_starts, y=y_starts, xend=x_ends, yend=y_ends)
    lines_df <- rbind(lines_df,square_lines_df)
  }
  
}

low_novel_allele_graph <- low_novel_allele_graph + annotate("segment", 
                                                            x = lines_df$x,
                                                            xend = lines_df$xend,
                                                            y = lines_df$y,
                                                            yend = lines_df$yend,
                                                            color = "gray50",
                                                            size=0.3)

plot(low_novel_allele_graph)


#############################################################################################################################
low_novel_appearance_df <- novel_appearance_df[novel_appearance_df$novel_allele %in% low_freq_allele_order,]

low_novel_appearance_df$carry <- unlist(lapply(1:nrow(low_novel_appearance_df), function(i){
  novel <- low_novel_appearance_df$novel_allele[i];
  subject <- low_novel_appearance_df$subject[i];
  subject <- sapply(strsplit(subject, "_", fixed = T), "[", 1);
  temp_row <- save_novel_geno[save_novel_geno$SUBJECT == subject & save_novel_geno$FULL_ALLELE_NAME == novel,];
  if (nrow(temp_row) == 0) {return("Not carry")};
  if_else(temp_row$HETERO, "Heterozygous", "Homozygous")
}))

low_snp_freq_plot <- ggplot(low_novel_appearance_df, aes(y=novel_allele, x=allele_freq, color = carry, shape=n_sequences)) + 
  geom_point() + theme_bw()+
  scale_y_discrete(limits = low_freq_allele_order) +
  scale_color_manual(values = c("skyblue", "blue", "grey")) +
  xlab('Allele frequency') + ylab('Polymorphism call') + 
  geom_hline(yintercept = 0.25, linetype="dashed", color = "gray50") +
  theme(axis.title = element_text(size = 16))
low_snp_freq_plot

low_snp_freq_plot <- low_snp_freq_plot + theme(axis.title = element_blank(), axis.text.y = element_blank(), legend.position = "none")
low_novel_allele_graph <- low_novel_allele_graph  + theme(legend.position = "none", axis.title = element_blank())


p1.common.y <- ggplot_gtable(ggplot_build(low_novel_allele_graph))
p2.common.y <- ggplot_gtable(ggplot_build(low_snp_freq_plot))
low_novel_x_width <- p1.common.y$widths

# copy the plot height from p1 to p2
p2.common.y$heights <- p1.common.y$heights

low_freq_novel_plot <- grid.arrange(p1.common.y, p2.common.y, ncol=2, widths=c(0.75,0.25))


################################################## GRAPH FOR HIGH ######################################################################################

# split the dataframe into 2 groups homozygous to the novel allele and heterozygous to the novel allele
homozygous_novel_df <- high_freq_novel_df
homozygous_novel_df$SUBJECTS_NUMBER <- homozygous_novel_df$HOMOZYGOUS
homozygous_novel_df$GROUP <- "Homozygous"
homozygous_novel_df <- homozygous_novel_df %>% select(GENE, FULL_ALLELE_NAME, SUBJECTS_NUMBER, GROUP)

heterozygous_novel_df <- high_freq_novel_df
heterozygous_novel_df$SUBJECTS_NUMBER <- heterozygous_novel_df$HETEROZYGOUS
heterozygous_novel_df$GROUP <- "Heterozygous"
heterozygous_novel_df <- heterozygous_novel_df %>% select(GENE, FULL_ALLELE_NAME, SUBJECTS_NUMBER, GROUP)

NOVEL_ALLELE_DF <- rbind(homozygous_novel_df, heterozygous_novel_df)
NOVEL_ALLELE_DF$FULL_ALLELE_NAME <- factor(NOVEL_ALLELE_DF$FULL_ALLELE_NAME, high_freq_allele_order)

colors <- unlist(lapply(high_freq_allele_order, function(allele) {
  if (novel_checks$both[novel_checks$ALLELE==allele]) {
    return("navy")
  } else if (novel_checks$padded[novel_checks$ALLELE==allele]) {
    return("green4")
  } else if ((novel_checks$unary[novel_checks$ALLELE==allele]) | (allele %in% extra_red_alleles)) {
    return("red")
  }
  return("black")
}))

high_novel_allele_graph <- ggplot(NOVEL_ALLELE_DF, aes(fill=GROUP, y=FULL_ALLELE_NAME, x=SUBJECTS_NUMBER)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("skyblue", "blue")) +
  scale_x_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(title="", reverse=TRUE)) +
  ylab("Undocumented pattern alleles")  + xlab("Number of individuals") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text.y= element_text(color = colors))

high_novel_allele_graph


lines_df <- c()
NOVEL_DIS_DF <- high_freq_novel_df
for (i in 1:nrow(NOVEL_DIS_DF)) {
  allele <- NOVEL_DIS_DF$FULL_ALLELE_NAME[i]
  y_ind <- which(high_freq_allele_order==allele)
  n_sides <- 3
  
  y_gaps <- c((1:(n_sides)/(n_sides+1))*0.84, rep(0, n_sides+1))
  y_starts <- ((y_ind-0.42)+rev(y_gaps))
  y_ends <- ((y_ind+0.42)-y_gaps)
  
  n_homozygous <- NOVEL_DIS_DF$HOMOZYGOUS[i]
  n_hetero <- NOVEL_DIS_DF$HETEROZYGOUS[i]
  
  if (NOVEL_DIS_DF$HOMOZYGOUS_LOW_K[i] > 0) {
    n_low_k <- NOVEL_DIS_DF$HOMOZYGOUS_LOW_K[i]
    x_gaps <- c((1:(n_sides)/(n_sides+1))*n_low_k, rep(0, n_sides+1))
    x_starts <- ((n_homozygous-n_low_k)+x_gaps)
    x_ends <- ((n_homozygous)-rev(x_gaps))
    square_lines_df <- data.frame(x=x_starts, y=y_starts, xend=x_ends, yend=y_ends)
    lines_df <- rbind(lines_df,square_lines_df)
  }
  
  if (NOVEL_DIS_DF$HETEROZYGOUS_LOW_K[i] > 0) {
    n_low_k <- NOVEL_DIS_DF$HETEROZYGOUS_LOW_K[i]
    x_gaps <- c((1:(n_sides)/(n_sides+1))*n_low_k, rep(0, n_sides+1))
    x_starts <- ((n_homozygous+n_hetero-n_low_k)+x_gaps)
    x_ends <- ((n_homozygous+n_hetero)-rev(x_gaps))
    square_lines_df <- data.frame(x=x_starts, y=y_starts, xend=x_ends, yend=y_ends)
    lines_df <- rbind(lines_df,square_lines_df)
  }
  
}

high_novel_allele_graph <- high_novel_allele_graph + annotate("segment", 
                                                              x = lines_df$x,
                                                              xend = lines_df$xend,
                                                              y = lines_df$y,
                                                              yend = lines_df$yend,
                                                              color = "gray50",
                                                              size=0.3)

plot(high_novel_allele_graph)


#############################################################################################################################
high_novel_appearance_df <- novel_appearance_df[novel_appearance_df$novel_allele %in% high_freq_allele_order,]

high_novel_appearance_df$carry <- unlist(lapply(1:nrow(high_novel_appearance_df), function(i){
  novel <- high_novel_appearance_df$novel_allele[i];
  subject <- high_novel_appearance_df$subject[i];
  subject <- sapply(strsplit(subject, "_", fixed = T), "[", 1);
  temp_row <- save_novel_geno[save_novel_geno$SUBJECT == subject & save_novel_geno$FULL_ALLELE_NAME == novel,];
  if (nrow(temp_row) == 0) {return("Not carry")};
  if_else(temp_row$HETERO, "Heterozygous", "Homozygous")
}))

high_snp_freq_plot <- ggplot(high_novel_appearance_df, aes(y=novel_allele, x=allele_freq, color = carry, shape=n_sequences)) + 
  geom_point() + theme_bw()+
  scale_y_discrete(limits = high_freq_allele_order) +
  scale_color_manual(values = c("skyblue", "blue", "grey")) +
  xlab('Allele frequency') + ylab('Polymorphism call') + 
  geom_hline(yintercept = 0.25, linetype="dashed", color = "gray50") +
  theme(axis.title = element_text(size = 16))
high_snp_freq_plot

high_snp_freq_plot <- high_snp_freq_plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
high_novel_allele_graph <- high_novel_allele_graph  + theme(axis.title.y = element_blank(), legend.position = "none")


p1.common.y <- ggplot_gtable(ggplot_build(high_novel_allele_graph))
p2.common.y <- ggplot_gtable(ggplot_build(high_snp_freq_plot))

p1.common.y$widths <- low_novel_x_width
# copy the plot height from p1 to p2
p2.common.y$heights <- p1.common.y$heights

high_freq_novel_plot <- grid.arrange(p1.common.y, p2.common.y, ncol=2, widths=c(0.75,0.25))

##############################################################################################################################################


# build legend
leg_df <- data.frame(x=rep(c(1:4), 4), y=unlist(lapply(1:4, rep, 4)), stringsAsFactors = F)
leg_df$carry <- ""
leg_texts <- c("", "lk < 3", "Undocumented allele not in genotype", "Heterozygous to the undocumented allele", "Homozygous to the undocumented allele")
leg_df$carry[leg_df$x == 1] <- leg_texts[-1]
leg_df$carry <- factor(leg_df$carry, levels = leg_texts)

leg_plot <- ggplot(leg_df, aes(x = x, y = y, fill = carry, label=carry)) +
  geom_tile() + 
  geom_text(nudge_x = 3.2) + 
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




g <- grid.arrange(grobs = list(low_freq_novel_plot, high_freq_novel_plot, leg_plot, shape_legend),
                  ylab = "Undocumented pattern alleles",
                  layout_matrix = rbind(c(1,1,1,NA),
                                        c(1,1,1,NA),
                                        c(1,1,1,NA),
                                        c(1,1,1,4),
                                        c(1,1,1,4),
                                        c(1,1,1,NA),
                                        c(1,1,1,3),
                                        c(1,1,1,3),
                                        c(1,1,1,NA),
                                        c(1,1,1,NA),
                                        c(2,2,2,NA),
                                        c(2,2,2,NA),
                                        c(2,2,2,NA)))

ggsave(paste0(figure_folder, "novel_DS4.pdf"),g, height = 8, width = 16)
ggsave(paste0(figure_folder, "novel_DS4.png"),g, height = 8, width = 16)
