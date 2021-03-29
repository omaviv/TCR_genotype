library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(gridExtra)
library(cowplot)
library(magick)
library(seqinr)

#####################################################################################################################
############################ Load required files ####################################################################
#####################################################################################################################
# the path to the repository
project_folder <- "tcr_genotype/"
required_files_folder <- paste0(project_folder, "figures/data/")
figure_folder <- paste0(project_folder, "figures/figure_1/")

# the distances between genes according to sequencing coverage
full_distances_df <- read.delim(paste0(required_files_folder, "TRBV_GENE_DISTANCES.tab"), sep = "\t", header=TRUE, row.names="GENES")
biomed_distances_df <- read.delim(paste0(required_files_folder, "BIOMED2_TRBV_GENE_DISTANCES.tab"), sep = "\t", header=TRUE, row.names="GENES")
adaptive_distances_df <- read.delim(paste0(required_files_folder, "Adaptive_TRBV_GENE_DISTANCES.tab"), sep = "\t", header=TRUE, row.names="GENES")

full_trbv <- read.fasta(paste0(required_files_folder, "TRBV_all_genes.fasta"),as.string = TRUE)
biomed_trbv <- read.fasta(paste0(required_files_folder, "BIOMED2_TRBV_all_genes.fasta"),as.string = TRUE)
adaptive_trbv <- read.fasta(paste0(required_files_folder, "Adaptive_TRBV_all_genes.fasta"),as.string = TRUE)

workflow <- ggdraw() + draw_image(paste0(figure_folder, "TCRB-workflow.png"))
vdj_seq <- ggdraw() + draw_image(paste0(figure_folder, "TCRB_Different_sequencing.png")) + theme(plot.margin = margin(0,10,0,50))

#####################################################################################################################
############################ General parameters #####################################################################
#####################################################################################################################

gray.colors(9, start = 0, end = 1, rev = FALSE)

brk <- c(0:5, 10, 25, 50)
labels <- c(as.character(0:4), "5-9", "10-24", "25-49", "50+")
labels_colors <- c('#FF0000FF','#FF3700FF','#FF4E00FF','#FF6F00FF','#FF8500FF','#FFA600FF','#FFD300FF','#FFFF70FF','#FFFFAFFF')
labels_colors <- gray.colors(9, start = 0, end = 0.9, rev = FALSE)

# find the gene order
distances_df <- as.matrix(full_distances_df)

rownames(distances_df) <- sub("TRB", "", rownames(distances_df))
colnames(distances_df) <- sub("TRB", "", colnames(distances_df))
colnames(distances_df) <- sub(".OR", "/OR", colnames(distances_df), fixed = T)
colnames(distances_df) <- gsub(".", "-", colnames(distances_df), fixed = T)

hcr <- hclust(dist(distances_df))
ddr <- as.dendrogram(hcr)
hcr$order

Rowv <- rowMeans(distances_df, na.rm = T)
ddr <- reorder(ddr, Rowv)

gene_order <- rownames(distances_df)[order.dendrogram(ddr)]

###############################################################################################################################
############################## FULL TRBV GENE DISTANCES #######################################################################
###############################################################################################################################

dt2 <- full_distances_df %>% 
  rownames_to_column() %>% 
  gather(colname, value, -rowname)

dt2$rowname <- sub("TRB", "", dt2$rowname)
dt2$colname <- sub("TRB", "", dt2$colname)
dt2$colname <- gsub(".OR", "/OR", dt2$colname, fixed = T)
dt2$colname <- gsub(".", "-", dt2$colname, fixed = T)

dt2$rowname <- factor(dt2$rowname, levels = gene_order)
dt2$colname <- factor(dt2$colname, levels = gene_order)

dt2$value[dt2$value >=50] <- 50
dt2$value[dt2$value >=25 & dt2$value <50] <- 25
dt2$value[dt2$value >=10 & dt2$value <25] <- 10
dt2$value[dt2$value >=5 & dt2$value <10] <- 5
dt2$value <- factor(dt2$value, levels = brk)

full_trbv_dis_graph <- ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()+ scale_fill_manual(values = labels_colors, breaks = brk, labels = labels, drop = FALSE) + 
  xlab("TRBV Genes") + ylab("TRBV Genes") +  labs(fill="Distance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8, angle = 90, vjust = 0.5), axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5), legend.position="bottom")

###############################################################################################################################
########################### BIOMED2 TRBV GENE DISTANCES #######################################################################
###############################################################################################################################

dt2 <- biomed_distances_df %>% 
  rownames_to_column() %>% 
  gather(colname, value, -rowname)

dt2$rowname <- sub("TRB", "", dt2$rowname)
dt2$colname <- sub("TRB", "", dt2$colname)
dt2$colname <- sub(".OR", "/OR", dt2$colname, fixed = T)
dt2$colname <- gsub(".", "-", dt2$colname, fixed = T)

dt2$rowname <- factor(dt2$rowname, levels = gene_order)
dt2$colname <- factor(dt2$colname, levels = gene_order)

dt2$value[dt2$value >=50] <- 50
dt2$value[dt2$value >=25 & dt2$value <50] <- 25
dt2$value[dt2$value >=10 & dt2$value <25] <- 10
dt2$value[dt2$value >=5 & dt2$value <10] <- 5
dt2$value <- factor(dt2$value, levels = brk)

biomed_trbv_dis_graph <- ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()+ scale_fill_manual(values = labels_colors, breaks = brk, labels = labels, drop = FALSE) + 
  xlab("TRBV Genes") + ylab("TRBV Genes") +  labs(fill="Distance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8, angle = 90, vjust = 0.5), axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5), legend.position="bottom")

###############################################################################################################################
########################### Adaptive TRBV GENE DISTANCES ######################################################################
###############################################################################################################################

dt2 <- adaptive_distances_df %>% 
  rownames_to_column() %>% 
  gather(colname, value, -rowname)

dt2$rowname <- sub("TRB", "", dt2$rowname)
dt2$colname <- sub("TRB", "", dt2$colname)
dt2$colname <- sub(".OR", "/OR", dt2$colname, fixed = T)
dt2$colname <- gsub(".", "-", dt2$colname, fixed = T)

dt2$rowname <- factor(dt2$rowname, levels = gene_order)
dt2$colname <- factor(dt2$colname, levels = gene_order)

dt2$value[dt2$value >=50] <- 50
dt2$value[dt2$value >=25 & dt2$value <50] <- 25
dt2$value[dt2$value >=10 & dt2$value <25] <- 10
dt2$value[dt2$value >=5 & dt2$value <10] <- 5
dt2$value <- factor(dt2$value, levels = brk)

adapt_trbv_dis_graph <- ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()+ scale_fill_manual(values = labels_colors, breaks = brk, labels = labels, drop = FALSE) + 
  xlab("TRBV Genes") + ylab("TRBV Genes") +  labs(fill="Distance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8, angle = 90, vjust = 0.5), axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5), legend.position="bottom")

#####################################################################################################################
############################ Unique sequences graph #################################################################
#####################################################################################################################

full_trbv <- as.character(names(full_trbv))
full_trbv <- sapply(base::strsplit(full_trbv, "*", fixed = T), "[", 1)
full_trbv <- as.data.frame(table(full_trbv))
names(full_trbv) <- c("Gene", "Allele_num")
full_trbv$sequencing <- "Full-length"

biomed_trbv <- as.character(names(biomed_trbv))
biomed_trbv <- sapply(base::strsplit(biomed_trbv, "*", fixed = T), "[", 1)
biomed_trbv <- as.data.frame(table(biomed_trbv))
names(biomed_trbv) <- c("Gene", "Allele_num")
biomed_trbv$sequencing <- "BIOMED-2"

adaptive_trbv <- as.character(names(adaptive_trbv))
adaptive_trbv <- sapply(base::strsplit(adaptive_trbv, "*", fixed = T), "[", 1)
adaptive_trbv <- as.data.frame(table(adaptive_trbv))
names(adaptive_trbv) <- c("Gene", "Allele_num")
adaptive_trbv$sequencing <- "Adaptive Biotechnologies"

gene_allele_num_df <- rbind(full_trbv, biomed_trbv, adaptive_trbv)
gene_allele_num_df$Gene <- gsub("TRB", "", gene_allele_num_df$Gene)
gene_allele_num_df$Gene <- factor(gene_allele_num_df$Gene, levels = gene_order)
gene_allele_num_df$sequencing <- factor(gene_allele_num_df$sequencing, levels = c("Full-length", "BIOMED-2", "Adaptive Biotechnologies"))

gene_allele_num_plot <- ggplot(data=gene_allele_num_df, aes(x=Gene, y=Allele_num, fill=sequencing)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_y_continuous(expand = c(0, 0)) +
  xlab("TRBV Genes") + ylab("Number of \n unique sequences") +  labs(fill="Sequencing approach") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=10, angle = 90, vjust = 0.5), axis.text.y = element_text(size=10))


#####################################################################################################################
############################ Collapse into figure ###################################################################
#####################################################################################################################
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(adapt_trbv_dis_graph)

full_trbv_dis_graph <- full_trbv_dis_graph + theme(legend.position="none")
biomed_trbv_dis_graph <- biomed_trbv_dis_graph + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank()) 
adapt_trbv_dis_graph <- adapt_trbv_dis_graph + theme(legend.position="none", axis.title.y = element_blank(), axis.text.y = element_blank())

require(grid)
a_label <- textGrob(
  label = "a.",
  gp = gpar(fontsize = 28), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

b_label <- textGrob(
  label = "b.",
  gp = gpar(fontsize = 28), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

c_label <- textGrob(
  label = "c.",
  gp = gpar(fontsize = 28), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

d_label <- textGrob(
  label = "d.",
  gp = gpar(fontsize = 28), 
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0)

workflow <- arrangeGrob(workflow, top = a_label)
vdj_seq <- arrangeGrob(vdj_seq, top = b_label)
distances_graph <- arrangeGrob(full_trbv_dis_graph, biomed_trbv_dis_graph, adapt_trbv_dis_graph, top = c_label, widths = c(1.15, 1, 1), ncol=3)
gene_allele_num_plot <- arrangeGrob(gene_allele_num_plot, top = d_label)

g <- grid.arrange(
  workflow, vdj_seq, distances_graph, mylegend,gene_allele_num_plot,
  heights = c(10, 10,11,1,5),
  layout_matrix = rbind(c(1),
                        c(2),
                        c(3),
                        c(4),
                        c(5))
)

ggsave(paste0(figure_folder, "Fig1.pdf"), plot = g, width = 19, height = 25)
ggsave(paste0(figure_folder, "Fig1.png"), plot = g, width = 19, height = 15)


# ggsave("TRBV_GENE_DISTANCES.pdf", plot = full_trbv_dis_graph, width = 8, height = 7)
# ggsave("TRBV_GENE_DISTANCES.png", plot = full_trbv_dis_graph, width = 8, height = 7)
# ggsave("BIOMED2_TRBV_GENE_DISTANCES.pdf", plot = biomed_trbv_dis_graph, width = 8, height = 7)
# ggsave("BIOMED2_TRBV_GENE_DISTANCES.png", plot = biomed_trbv_dis_graph, width = 8, height = 7)
# ggsave("Adaptive_TRBV_GENE_DISTANCES.pdf", plot = adapt_trbv_dis_graph, width = 8, height = 7)
# ggsave("Adaptive_TRBV_GENE_DISTANCES.png", plot = adapt_trbv_dis_graph, width = 8, height = 7)
