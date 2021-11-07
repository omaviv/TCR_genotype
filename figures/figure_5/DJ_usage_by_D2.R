library(ggplot2)
library(gplots)
library(dplyr)
library(data.table)
library(fastmatch)
library(reshape2)
library(ggsignif)
library(gridExtra)

##################################################################################################################################################
############################ Load required files #################################################################################################
##################################################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_5/")
references_folder <- paste0(project_folder, "pipeline/fasta_references/")

# Load TRB locus and pseudo genes
load(paste0(project_folder, "pipeline/sysdata.rda"))

pairing_df <- read.delim(paste0(required_files_folder, "DS4_all_dj_pairs.tab"), sep = "\t", stringsAsFactors = F)
genotypes <- read.delim(paste0(required_files_folder, "DS4_All_Genotypes.tab"), sep = "\t", stringsAsFactors = F)
trbj_usage <- read.delim(paste0(required_files_folder, "Adaptive_TRBJ_gene_usage.tab"), header = T, sep = "\t", stringsAsFactors = F)

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

x_space <- 0.25
asterisk_size <- 3

TRBD2_geno <- genotypes[genotypes$GENE=="TRBD2",]
TRBD2_geno$GENOTYPED_ALLELES[grepl(",", TRBD2_geno$GENOTYPED_ALLELES)] <- "01,02"
TRBD2_geno$d2_geno <- TRBD2_geno$GENOTYPED_ALLELES
TRBD2_geno <- TRBD2_geno %>% select(SUBJECT, d2_geno)
names(TRBD2_geno) <- tolower(names(TRBD2_geno))
pairing_df <- merge(pairing_df, TRBD2_geno)

pairing_df <- pairing_df[pairing_df$d_len >= 8 & pairing_df$productive,]

trbj_usage <- trbj_usage[!grepl("P|,", trbj_usage$j_gene) & trbj_usage$productive,]

trbj_usage <- trbj_usage %>% group_by(subject, j_gene) %>% dplyr::summarise(N=sum(N))
trbj_usage <- merge(trbj_usage, TRBD2_geno)

TRBJ_GENES <- GENE.loc[["TRB"]]
TRBJ_GENES <- TRBJ_GENES[grepl("J[1-2]-[1-7]", TRBJ_GENES)]

gene_colors <- ifelse(grepl("J2-", TRBJ_GENES), "blue", "red")
names(gene_colors) <- TRBJ_GENES

##################################################################################################################################################
################### TRBD & TRBJ usage by TRBD2 genotype ##########################################################################################
##################################################################################################################################################

d_usage_df <- pairing_df[!grepl(",", pairing_df$d_gene),]
total_d <- d_usage_df %>% group_by(subject) %>% dplyr::summarise(total=sum(number_of_comb))
d_usage_df <- d_usage_df %>% group_by(subject, d_gene, d2_geno) %>% dplyr::summarise(N=sum(number_of_comb))
d_usage_df <- merge(d_usage_df, total_d)
d_usage_df$usage <- d_usage_df$N / d_usage_df$total

d2_usage_graph <- ggplot(d_usage_df, aes(x=gsub("TRB", "", d_gene), y=usage, fill = d2_geno)) + 
  geom_boxplot(outlier.size = 0.4) +  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(d_usage_df$usage) *1.1)) +
  ylab("Usage") + xlab("") +
  labs(fill="TRBD2 genotype") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

temp <- d_usage_df[d_usage_df$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$d_gene), function(x){wilcox.test(usage~d2_geno, x)$p.value})

labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
d2_usage_graph <- d2_usage_graph + geom_signif(y_position = max(d_usage_df$usage), xmin = 1:2 - x_space, xmax = 1:2 + x_space, annotations = labels, textsize = asterisk_size)

d2_usage_graph

trbj_total_usage <- trbj_usage %>% group_by(subject) %>% dplyr::summarise(total = sum(N))
trbj_usage <- merge(trbj_usage, trbj_total_usage)
trbj_usage$usage <- trbj_usage$N / trbj_usage$total

trbj_usage_graph <- ggplot(trbj_usage, aes(x=gsub("TRB", "", j_gene), y=usage, fill = d2_geno)) + 
  geom_boxplot(outlier.size = 0.4) +  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(trbj_usage$usage) *1.1)) +
  ylab("Usage") + xlab("TRBJ gene") +
  guides(fill=guide_legend("TRBD2 genotype", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.5,0.9), axis.text.x=element_text(color=gene_colors, size=10),
        legend.background = element_blank(), legend.box.background = element_rect(colour = "black", fill = "grey90"))



# y.values <- c(0.32, 0.325, 0.35)
y.values <- sapply(split(trbj_usage, trbj_usage$j_gene), function(x){max(x$usage)}) + 0.01

temp <- trbj_usage[trbj_usage$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_gene), function(x){wilcox.test(usage~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
trbj_usage_graph <- trbj_usage_graph + geom_signif(y_position = y.values, xmin = 1:13 - x_space, xmax = 1:13 + x_space, annotations = labels, textsize = asterisk_size)
top_legend <- g_legend(trbj_usage_graph)


a <- d2_usage_graph  + theme(legend.position="none", axis.title.x = element_blank())
b <- trbj_usage_graph + theme(legend.position="none", axis.title.y = element_blank(), axis.title.x = element_blank())
adapt_d2_geno_leg<-g_legend(d2_usage_graph)


adaptive_trbj_usage <- arrangeGrob(a, b, adapt_d2_geno_leg, bottom = "Gene", widths = c(0.25,0.6,0.15), ncol=3)

# ggsave(paste0(figure_folder, "DS4_TRBJ_USAGE.png"), adaptive_trbj_usage, height = 6)
# ggsave(paste0(figure_folder, "DS4_TRBJ_USAGE.pdf"), adaptive_trbj_usage, width = 20, height = 6)

##################################################################################################################################################
################### Total TRBJ1 & TRBJ2 usage by TRBD2 genotype #############################################################################
##################################################################################################################################################
trbj_family_usage <- trbj_usage
trbj_family_usage$j_family <- substr(trbj_family_usage$j_gene, 1, 5)
trbj_family_usage <- trbj_family_usage %>% group_by(subject, j_family, d2_geno) %>% dplyr::summarise(total_family = sum(N))
trbj_total_usage <- trbj_family_usage %>% group_by(subject) %>% dplyr::summarise(total = sum(total_family))
trbj_family_usage <- merge(trbj_family_usage, trbj_total_usage)
trbj_family_usage$total_usage <- trbj_family_usage$total_family / trbj_family_usage$total

# calculate p values
temp <- trbj_family_usage[trbj_family_usage$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_family), function(x){wilcox.test(total_usage~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

trbj_family_usage_graph <- ggplot(trbj_family_usage, aes(x=gsub("TRB", "", j_family), y=total_usage, fill = d2_geno)) + 
  geom_boxplot(outlier.size = 0.4) +  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(trbj_family_usage$total_usage) *1.3)) +
  ylab("Usage") + xlab("TRBJ family") +
  guides(fill=guide_legend("TRBD2 genotype", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        legend.position=c(0.5,0.9))
# legend.position=c(0.5,0.9), axis.text.x=element_text(color=gene_colors))

# y.values <- c(0.32, 0.325, 0.35)
y.values <- sapply(split(trbj_family_usage, trbj_family_usage$j_family), function(x){max(x$total_usage)}) 
y.values <- lapply(y.values, function(x){max(x*1.05, x+0.02)})
y.values <- unlist(y.values)

labels <- symnum(p.values[grepl("J1", names(p.values))], corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
trbj_family_usage_graph <- trbj_family_usage_graph + geom_signif(y_position = y.values, xmin = 1:length(y.values) - x_space, xmax = 1:length(y.values) + x_space, annotations = labels, textsize = asterisk_size)
trbj_family_usage_graph

# ggsave(paste0(figure_folder, "DS4_TRBJ_family_usage.png"), trbj_family_usage_graph, height = 6)
# ggsave(paste0(figure_folder, "DS4_TRBJ_family_usage.pdf"), trbj_family_usage_graph)


##################################################################################################################################################
################### normalized TRBJ1 & TRBJ2 usage by TRBD2 genotype #############################################################################
##################################################################################################################################################


trbj_usage$j_family <- substr(trbj_usage$j_gene, 1, 5)
trbj_family_total <- trbj_usage %>% group_by(subject, j_family) %>% dplyr::summarise(total_family = sum(N))
trbj_usage <- merge(trbj_usage, trbj_family_total)
trbj_usage$norm_usage <- trbj_usage$N / trbj_usage$total_family

# calculate p values
temp <- trbj_usage[trbj_usage$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_gene), function(x){wilcox.test(norm_usage~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

trbj1_usage <- trbj_usage[grepl("J1", trbj_usage$j_gene),]
trbj1_usage_graph <- ggplot(trbj1_usage, aes(x=gsub("TRB", "", j_gene), y=norm_usage, fill = d2_geno)) + 
  geom_boxplot(outlier.size = 0.4) +  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(trbj_usage$norm_usage) *1.1)) +
  ylab("P(TRBJ1-N)") + xlab("TRBJ1 gene") +
  guides(fill=guide_legend("TRBD2 genotype", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        legend.position=c(0.5,0.9))
# legend.position=c(0.5,0.9), axis.text.x=element_text(color=gene_colors))

# y.values <- c(0.32, 0.325, 0.35)
y.values <- sapply(split(trbj1_usage, trbj1_usage$j_gene), function(x){max(x$norm_usage)}) 
y.values <- lapply(y.values, function(x){max(x*1.05, x+0.02)})
y.values <- unlist(y.values)

labels <- symnum(p.values[grepl("J1", names(p.values))], corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
trbj1_usage_graph <- trbj1_usage_graph + geom_signif(y_position = y.values, xmin = 1:length(y.values) - x_space, xmax = 1:length(y.values) + x_space, annotations = labels, textsize = asterisk_size)


trbj2_usage <- trbj_usage[grepl("J2", trbj_usage$j_gene),]
trbj2_usage_graph <- ggplot(trbj2_usage, aes(x=gsub("TRB", "", j_gene), y=norm_usage, fill = d2_geno)) + 
  geom_boxplot(outlier.size = 0.4) +  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(trbj_usage$norm_usage) *1.1)) +
  ylab("P(TRBJ2-N)") + xlab("TRBJ2 gene") +
  guides(fill=guide_legend("TRBD2 genotype", nrow = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        legend.position=c(0.5,0.9))
# legend.position=c(0.5,0.9), axis.text.x=element_text(color=gene_colors))

# y.values <- c(0.32, 0.325, 0.35)
y.values <- sapply(split(trbj2_usage, trbj2_usage$j_gene), function(x){max(x$norm_usage)}) 
y.values <- lapply(y.values, function(x){max(x*1.05, x+0.02)})
y.values <- unlist(y.values)

labels <- symnum(p.values[grepl("J2", names(p.values))], corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
trbj2_usage_graph <- trbj2_usage_graph + geom_signif(y_position = y.values, xmin = 1:length(y.values) - x_space, xmax = 1:length(y.values) + x_space, annotations = labels, textsize = asterisk_size)

a <- trbj1_usage_graph  + theme(legend.position="none", axis.title.x = element_blank())
b <- trbj2_usage_graph + theme(legend.position="none", axis.text.y = element_blank(), axis.title.x = element_blank())

trbj_norm_usage_graph <- arrangeGrob(a, b, adapt_d2_geno_leg, bottom = "J gene", widths = c(0.5, 0.5,0.15), ncol=3)
plot(trbj_norm_usage_graph)

# ggsave(paste0(figure_folder, "DS4_norm_TRBJ_USAGE.png"), trbj_norm_usage_graph, height = 6)
# ggsave(paste0(figure_folder, "DS4_norm_TRBJ_USAGE.pdf"), trbj_norm_usage_graph)


##################################################################################################################################################
################### TRBD2-TRBJ2 pairing ##########################################################################################################
##################################################################################################################################################

d2_j2_pairing <- pairing_df[grepl("J2", pairing_df$j_gene) & pairing_df$d_gene=="TRBD2",] 
d2_j2_pairing <- d2_j2_pairing %>% group_by(subject, d2_geno, d_gene, j_gene) %>% dplyr::summarise(number_of_comb=sum(number_of_comb))

d2_total <- d2_j2_pairing %>% group_by(subject, d_gene) %>% dplyr::summarise(total_d2=sum(number_of_comb))

d2_j2_pairing <- merge(d2_j2_pairing, d2_total, by=c("subject", "d_gene"))
d2_j2_pairing$frac <- d2_j2_pairing$number_of_comb / d2_j2_pairing$total_d2
d2_j2_pairing$j_gene <- gsub("TRB", "", d2_j2_pairing$j_gene)
d2_j2_pairing <- d2_j2_pairing[!grepl(",|P", d2_j2_pairing$j_gene),]

d2_j2_pairing_graph <- ggplot(d2_j2_pairing, aes(x=j_gene, y=frac)) + 
  geom_boxplot(aes(fill=d2_geno)) + scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(d2_j2_pairing$frac) *1.1)) +
  xlab("TRBJ2 gene") + ylab("P(TRBJ2-N|TRBD2)") + 
  guides(fill=guide_legend("TRBD2 genotype")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

d2_j2_pairing_graph
y.values <- sapply(split(d2_j2_pairing, d2_j2_pairing$j_gene), function(x){max(x$frac)}) + max(d2_j2_pairing$frac)/20

temp <- d2_j2_pairing[d2_j2_pairing$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_gene), function(x){wilcox.test(frac~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
d2_j2_pairing_graph <- d2_j2_pairing_graph + geom_signif(y_position = y.values, xmin = 1:7 - x_space, xmax = 1:7 + x_space, annotations = labels, textsize = asterisk_size)

# ggsave(paste0(figure_folder, "DS4_TRBD2_J2_pairing.pdf"), d2_j2_pairing_graph)
# ggsave(paste0(figure_folder, "DS4_TRBD2_J2_pairing.png"), d2_j2_pairing_graph)

##################################################################################################################################################
################### TRBD1-TRBJ familty pairing ###################################################################################################
##################################################################################################################################################

d1_j_fam_pairing <- pairing_df[pairing_df$d_gene=="TRBD1",] 
d1_j_fam_pairing$j_family <- "TRBJ1"
d1_j_fam_pairing$j_family[grepl("J2", d1_j_fam_pairing$j_gene)] <- "TRBJ2"
d1_j_fam_pairing <- d1_j_fam_pairing %>% group_by(subject, d2_geno, d_gene, j_family) %>% dplyr::summarise(number_of_comb=sum(number_of_comb))

d1_total <- d1_j_fam_pairing %>% group_by(subject, d_gene) %>% dplyr::summarise(total_d1=sum(number_of_comb))

d1_j_fam_pairing <- merge(d1_j_fam_pairing, d1_total, by=c("subject", "d_gene"))
d1_j_fam_pairing$frac <- d1_j_fam_pairing$number_of_comb / d1_j_fam_pairing$total_d1
d1_j_fam_pairing <- d1_j_fam_pairing[!grepl(",", d1_j_fam_pairing$j_family),]

d1_j_fam_pairing_graph <- ggplot(d1_j_fam_pairing, aes(x=gsub("TRB", "", j_family), y=frac)) + 
  geom_boxplot(aes(fill=d2_geno)) + scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(d1_j_fam_pairing$frac) *1.1)) +
  xlab("TRBJ family") + ylab("P(TRBJ|TRBD1)") + 
  guides(fill=guide_legend("TRBD2 genotype")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

y.values <- sapply(split(d1_j_fam_pairing, d1_j_fam_pairing$j_family), function(x){max(x$frac)}) + max(d1_j_fam_pairing$frac)/20

temp <- d1_j_fam_pairing[d1_j_fam_pairing$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_family), function(x){wilcox.test(frac~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
d1_j_fam_pairing_graph <- d1_j_fam_pairing_graph + geom_signif(y_position = y.values, xmin = 1:2 - x_space, xmax = 1:2 + x_space, annotations = labels, textsize = asterisk_size)

# ggsave(paste0(figure_folder, "DS4_TRBD1_J_family_pairing.pdf"), d1_j_fam_pairing_graph)
# ggsave(paste0(figure_folder, "DS4_TRBD1_J_family_pairing.png"), d1_j_fam_pairing_graph)

##################################################################################################################################################
################### TRBD1-TRBJ pairing ##########################################################################################################
##################################################################################################################################################

d1_j_pairing <- pairing_df[pairing_df$d_gene=="TRBD1",] 
d1_j_pairing <- d1_j_pairing %>% group_by(subject, d2_geno, d_gene, j_gene) %>% dplyr::summarise(number_of_comb=sum(number_of_comb))

d1_total <- d1_j_pairing %>% group_by(subject, d_gene) %>% dplyr::summarise(total_d1=sum(number_of_comb))

d1_j_pairing <- merge(d1_j_pairing, d1_total, by=c("subject", "d_gene"))
d1_j_pairing$frac <- d1_j_pairing$number_of_comb / d1_j_pairing$total_d1
d1_j_pairing$j_gene <- gsub("TRB", "", d1_j_pairing$j_gene)
d1_j_pairing <- d1_j_pairing[!grepl(",|P", d1_j_pairing$j_gene),]

d1_j_pairing_graph <- ggplot(d1_j_pairing, aes(x=j_gene, y=frac)) + 
  geom_boxplot(aes(fill=d2_geno)) + scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(d1_j_pairing$frac) *1.1)) +
  guides(fill=guide_legend("TRBD2 genotype")) +
  xlab("TRBJ gene") + ylab("P(TRBJ1/2-N|TRBD1)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        axis.text.x=element_text(size=10, color=gene_colors))

y.values <- sapply(split(d1_j_pairing, d1_j_pairing$j_gene), function(x){max(x$frac)}) + max(d1_j_pairing$frac)/20

temp <- d1_j_pairing[d1_j_pairing$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_gene), function(x){wilcox.test(frac~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
d1_j_pairing_graph <- d1_j_pairing_graph + geom_signif(y_position = y.values, xmin = 1:13 - x_space, xmax = 1:13 + x_space, annotations = labels, textsize = asterisk_size)

# ggsave(paste0(figure_folder, "DS4_TRBD1_J_genes_pairing.pdf"), d1_j_pairing_graph)
# ggsave(paste0(figure_folder, "DS4_TRBD1_J_genes_pairing.png"), d1_j_pairing_graph)

##################################################################################################################################################
################### TRBD1-TRBJ2 pairing ##########################################################################################################
##################################################################################################################################################

d1_j2_pairing <- pairing_df[grepl("J2", pairing_df$j_gene) & pairing_df$d_gene=="TRBD1",] 
d1_j2_pairing <- d1_j2_pairing %>% group_by(subject, d2_geno, d_gene, j_gene) %>% dplyr::summarise(number_of_comb=sum(number_of_comb))

d1_total <- d1_j2_pairing %>% group_by(subject, d_gene) %>% dplyr::summarise(total_d1=sum(number_of_comb))

d1_j2_pairing <- merge(d1_j2_pairing, d1_total, by=c("subject", "d_gene"))
d1_j2_pairing$frac <- d1_j2_pairing$number_of_comb / d1_j2_pairing$total_d1
d1_j2_pairing$j_gene <- gsub("TRB", "", d1_j2_pairing$j_gene)
d1_j2_pairing <- d1_j2_pairing[!grepl(",|P", d1_j2_pairing$j_gene),]

d1_j2_pairing_graph <- ggplot(d1_j2_pairing[!grepl(",", d1_j2_pairing$j_gene),], aes(x=j_gene, y=frac)) + 
  geom_boxplot(aes(fill=d2_geno)) + scale_fill_brewer(palette = "Dark2") +
  xlab("TRBJ2 gene") + ylab("P(TRBJ2-N|TRBD1)") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(d1_j2_pairing$frac) *1.1)) +
  guides(fill=guide_legend("TRBD2 genotype")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10),
        axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

d1_j2_pairing_graph
y.values <- sapply(split(d1_j2_pairing, d1_j2_pairing$j_gene), function(x){max(x$frac)}) + max(d1_j2_pairing$frac)/20

temp <- d1_j2_pairing[d1_j2_pairing$d2_geno %in% c("01", "02"),] 
p.values <- sapply(split(temp, temp$j_gene), function(x){wilcox.test(frac~d2_geno, x)$p.value})
p.values <- p.adjust(p.values, method = "bonferroni")

labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*", "n.s."))
d1_j2_pairing_graph <- d1_j2_pairing_graph + geom_signif(y_position = y.values, xmin = 1:7 - x_space, xmax = 1:7 + x_space, annotations = labels, textsize = asterisk_size)

# ggsave(paste0(figure_folder, "DS4_TRBD1_J2_pairing.pdf"), d1_j2_pairing_graph)
# ggsave(paste0(figure_folder, "DS4_TRBD1_J2_pairing.png"), d1_j2_pairing_graph)

##################################################################################################################################################
################### Merge all 6 graphs into 1 figure #############################################################################################
##################################################################################################################################################

a_label <- textGrob(
  label = "A.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

b_label <- textGrob(
  label = "B.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

c_label <- textGrob(
  label = "C.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

d_label <- textGrob(
  label = "D.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

e_label <- textGrob(
  label = "E.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

f_label <- textGrob(
  label = "F.",
  gp = gpar(fontsize = 20), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

d2_usage_graph <- d2_usage_graph  + theme(legend.position="none", axis.title.x = element_blank())
trbj_usage_graph <- trbj_usage_graph + theme(legend.position="none", axis.title.y = element_blank(), axis.title.x = element_blank())
trbdj_usage_graph <- arrangeGrob(d2_usage_graph, trbj_usage_graph, bottom = textGrob(label = "TRB gene", gp = gpar(fontsize = 16)), widths = c(0.25,0.75), 
                                 ncol=2, top = a_label)

trbj_family_usage_graph <- trbj_family_usage_graph + theme(legend.position="none")
trbj_family_usage_graph <- arrangeGrob(trbj_family_usage_graph, top = b_label)

trbj1_usage_graph <- trbj1_usage_graph + theme(legend.position="none", axis.title.x = element_blank())
trbj2_usage_graph <- trbj2_usage_graph + theme(legend.position="none", axis.text.y = element_blank(), axis.title.x = element_blank())
trbj_norm_usage_graph <- arrangeGrob(trbj1_usage_graph, trbj2_usage_graph, bottom = textGrob(label = "TRBJ gene", gp = gpar(fontsize = 16)), widths = c(0.63, 0.67), ncol=2, top = c_label)


d2_j2_pairing_graph <- d2_j2_pairing_graph + theme(legend.position="none")
d2_j2_pairing_graph <- arrangeGrob(d2_j2_pairing_graph, top = d_label)


d1_j_pairing_graph <- d1_j_pairing_graph + theme(legend.position="none")
d1_j_pairing_graph <- arrangeGrob(d1_j_pairing_graph, top = e_label)

d1_j2_pairing_graph <- d1_j2_pairing_graph + theme(legend.position="none")
d1_j2_pairing_graph <- arrangeGrob(d1_j2_pairing_graph, top = f_label)


g <- grid.arrange(
  grobs = list(trbdj_usage_graph, trbj_family_usage_graph, trbj_norm_usage_graph, d2_j2_pairing_graph,
               d1_j_pairing_graph, d1_j2_pairing_graph),
  heights = c(0.3,0.3,0.3, 0.03),
  nrow = 4,
  bottom = top_legend,
  padding = unit(2, "line"),
  layout_matrix = rbind(
    c(1, 4),
    c(2, 5),
    c(3, 6), 
    c())
)

ggsave(paste0(figure_folder, "Fig5.pdf"), g, height = 11, width = 13)
# ggsave(paste0(figure_folder, "D_J_usage_prod_len_8.png"), g, height = 10, width = 12)

# d1_j_pairing_med <- d1_j_pairing %>% group_by(d2_geno, j_gene) %>% summarise(median_frac = median(frac))
# d1_j_pairing_med <- d1_j_pairing_med[!grepl(",", d1_j_pairing_med$d2_geno),]
# d1_j_pairing_med <- d1_j_pairing_med[grepl("J2", d1_j_pairing_med$j_gene),]
# 
# for (j in unique(d1_j_pairing_med$j_gene)) {
#   temp <- d1_j_pairing_med[d1_j_pairing_med$j_gene==j,]
#   a <- temp$median_frac[temp$d2_geno == "02"] / temp$median_frac[temp$d2_geno == "01"]
#   print(paste(j, a))
# }
# 
# d1_j2_pairing_med <- d1_j2_pairing %>% group_by(d2_geno, j_gene) %>% summarise(median_frac = median(frac))
# d1_j2_pairing_med <- d1_j2_pairing_med[!grepl(",", d1_j2_pairing_med$d2_geno),]
# 
# for (j in unique(d1_j2_pairing_med$j_gene)) {
#   temp <- d1_j2_pairing_med[d1_j2_pairing_med$j_gene==j,]
#   a <- temp$median_frac[temp$d2_geno == "02"] / temp$median_frac[temp$d2_geno == "01"]
#   print(paste(j, a))
# }
