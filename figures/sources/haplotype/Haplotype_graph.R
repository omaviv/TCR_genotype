# HaplotypeR graphic functions -----------------------------------------------------

#' @include rabhit.R
#' @include internal_functions.R
NULL

#' Graphical output of an inferred haplotype
#'
#' The \code{plotHaplotype} functions visualizes an inferred haplotype.
#'
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    html_output          if TRUE, a html5 interactive graph is outputed. Defualt is FALSE.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    text_size            the size of graph labels. Default is 14 (pts).
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    plotYaxis            if TRUE, Y axis labels (gene names) are plotted on the middle and right plots. Default is TRUE.
#' @param    chain                the Ig chain: IGH,IGK,IGL. Default is IGH.
#' @param    dir                  The output folder for saving the haplotype map for multiple individuals.
#'
#' @return
#'
#' A haplotype map visualization. If more than one subject is visualized, a pdf is created. If html_output is TRUE, a folder named html_output is created with individual graphs.
#'
#' @details
#'
#' A \code{data.frame} in a haplotype format created by \code{createFullHaplotype} function.
#'
#' @examples
#'
#' # Selecting a single individual from the haplotype samples data
#' haplo_db = samplesHaplotype[samplesHaplotype$SUBJECT=='I5', ]
#'
#' # plot haplotype
#' plotHaplotype(haplo_db)
#'
#' @export
plotHaplotype <- function(hap_table, html_output = FALSE, gene_sort = c("name", "position"), text_size = 14, removeIGH = TRUE, plotYaxis = TRUE, chain = c("IGH","IGK", "IGL", "TRB"), dir) {
  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)
  
  if (missing(gene_sort)) {
    gene_sort = "position"
  }
  gene_sort <- match.arg(gene_sort)
  
  id_GENE_col <- which(names(hap_table)=="GENE")
  hapBy_col_id <- c(id_GENE_col+1,id_GENE_col+2)
  hapBy_cols = names(hap_table)[c(id_GENE_col+1,id_GENE_col+2)]
  
  hapBy_alleles = gsub("_", "*", hapBy_cols)
  
  if (!("SUBJECT" %in% names(hap_table))) {
    hap_table$SUBJECT <- rep("S1", nrow(hap_table))
  }
  
  plot_list <- c()
  for (sample_name in unique(hap_table$SUBJECT)) {
    
    
    GENE.loc.tmp <- GENE.loc[[chain]]
    
    haplo.db <- parseHapTab(hap_table[hap_table$SUBJECT == sample_name, ], chain = chain, sample_name = sample_name,
                            hapBy_cols = hapBy_cols, hapBy_alleles = hapBy_alleles)
    geno.df <- sortDFByGene(haplo.db$geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    kval.df <- sortDFByGene(haplo.db$kval.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    count.df <- sortDFByGene(haplo.db$count.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    
    ########################################################################################################
    
    ### Prepare All panels
    geno.df$ALLELE_TEXT <- geno.df$ALLELES
    count.df$ALLELE_TEXT <- count.df$ALLELES
    
    
    if (length(grep("[0-9][0-9]_[0-9][0-9]$", geno.df$ALLELES)) != 0) {
      geno.df <- as.data.frame(geno.df %>% group_by(.data$hapBy, .data$GENE) %>% mutate(n = dplyr::n()))
      geno.df$freq <- ifelse(geno.df$n == 2, 0.5, ifelse(geno.df$n != 1, 0.25, 1))
      non_reliable_alleles_text <- nonReliableAllelesText_V2(non_reliable_alleles_text = geno.df[grep("[0-9][0-9]_[0-9][0-9]$", geno.df$ALLELES), ])
    } else {
      non_reliable_alleles_text <- c()
    }
    
    geno.df$ALLELE_TEXT <- sapply(1:nrow(geno.df),function(i){
      if(grepl("[0-9][0-9]_[0-9][0-9]$", geno.df$ALLELES[i])){
        unique(non_reliable_alleles_text$text_bottom[geno.df$GENE[i]==non_reliable_alleles_text$GENE & grepl(geno.df$ALLELES[i],non_reliable_alleles_text$text_bottom)])
      }else{
        geno.df$ALLELES[i]
      }
    })
    
    count.df$ALLELE_TEXT <- sapply(1:nrow(count.df),function(i){
      if(grepl("[0-9][0-9]_[0-9][0-9]$", count.df$ALLELES[i])){
        unique(non_reliable_alleles_text$text_bottom[count.df$GENE[i]==non_reliable_alleles_text$GENE & grepl(count.df$ALLELES[i],non_reliable_alleles_text$text_bottom)])
      }else{
        count.df$ALLELES[i]
      }
    })
    
    ## Middle panel
    geno.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]$", geno.df$ALLELES)] <- "NRA"
    
    allele_palette <- alleleHapPalette(geno.df$ALLELES)
    AlleleCol <- allele_palette$AlleleCol
    transper <- allele_palette$transper
    
    geno.df$ALLELES <- factor(geno.df$ALLELES, levels = AlleleCol)
    geno.df$text <- paste0("Gene: ", geno.df$GENE, '<br />',"Allele: ", geno.df$ALLELE_TEXT)
    options(warn = -1)
    p = ggplot() +
      geom_bar(data = geno.df, mapping = aes_string(x = "GENE", fill = "ALLELES", text = "text"),
               position = "fill", width = 0.9, na.rm = T) + coord_flip() +
      xlab("") + ylab("") + facet_grid(paste0(".~", "hapBy"), switch = "x") + theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
                                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = text_size), strip.background = element_blank(),
                                                                                                 strip.text = element_text(face = "bold"), axis.text = element_text(colour = "black"), panel.spacing = unit(0, "cm"), strip.switch.pad.grid = unit(0,
                                                                                                                                                                                                                                                   "cm"), plot.margin = unit(c(0.25, 0, 0.2, 0), "cm"), legend.key = element_rect("#DCDCDC")) + scale_fill_manual(values = alpha(names(AlleleCol), transper), name = "Alleles", drop = F)
    
    if (!plotYaxis) {
      p = p + theme(axis.text.y = element_blank())
    }
    
    ## Right panel plot K values
    kval.df$text <- paste0("Gene: ", kval.df$GENE, '<br />',"lK: ", round(as.numeric(kval.df$K),4), '<br />',"lK group: ",kval.df$K_GROUPED)
    pk <- ggplot(kval.df, aes_string(x = "GENE", fill = "K_GROUPED", text = "text")) + theme_bw() + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
                                                                                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = text_size), strip.background = element_blank(),
                                                                                                          strip.text = element_text(face = "bold"), panel.spacing = unit(0, "cm"), strip.switch.pad.grid = unit(0, "cm"), plot.margin = unit(c(0.25, 0,
                                                                                                                                                                                                                                               0.2, 0), "cm"), legend.key = element_rect("#DCDCDC")) + geom_bar(position = "fill", width = 0.7, na.rm = T) + coord_flip() + xlab("") + ylab("") + facet_grid(paste0(".~", "hapBy"),
                                                                                                                                                                                                                                                                                                                                                                                                             switch = "x")
    
    count.df$ALLELES <- as.character(count.df$ALLELES)
    count.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]$", count.df$ALLELES)] <- "NRA"
    count.df$border <- factor(ifelse(count.df$ALLELES == "NRA", "black", "white"), levels = c("black", "white"))
    count.df$ALLELES <- factor(count.df$ALLELES, levels = AlleleCol)
    count.df$text <- paste0("Gene: ", count.df$GENE, '<br />',"Allele: ", count.df$ALLELE_TEXT,
                            '<br />',"Count: ", round(count.df$COUNT,4),
                            '<br />',"log<sub>10</sub>(Count+1):", abs(round(count.df$COUNT2,4)))
    ## Left panel
    p2 <- ggplot(count.df, aes_string(x = "GENE", y = "COUNT2", fill = "ALLELES",
                                      text = "text")) +
      geom_bar(stat = "identity", position = "Dodge", width = 0.9,
               na.rm = T, aes_string(colour = "border")) + coord_flip() + cowplot::background_grid(minor = "none") +
      scale_fill_manual(values = alpha(names(AlleleCol),
                                       transper), name = "ALLELES", drop = F) + scale_color_manual(values = alpha(c("black", "white"), c(0.5, 0)), drop = F) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"), axis.text = element_text(colour = "black"),
            text = element_text(size = text_size), plot.margin = unit(c(0.25,
                                                                        0, -0.05, 0), "cm"), panel.background = element_blank(), legend.key = element_rect("#DCDCDC")) +
      scale_y_continuous(breaks = seq(-3, 3, by = 1), labels = c(3:0, 1:3)) + ylab(expression("log"[10] *
                                                                                                "(Count+1)")) + xlab("Gene") + geom_hline(yintercept = c(0), linetype = "dotted")
    
    short_reads = F
    if (is.data.frame(non_reliable_alleles_text)) {
      p <- p + geom_text(data = non_reliable_alleles_text, aes_string(label = "text", x = "GENE", y = "pos"),
                         angle = 0, size = non_reliable_alleles_text$size)
      non_reliable_alleles_text$COUNT2 <- count.df$COUNT2[count.df$ALLELES == "NRA"]
      non_reliable_alleles_text <- non_reliable_alleles_text[non_reliable_alleles_text$COUNT2 != 0, ]
      non_reliable_alleles_text$hjust <- ifelse(non_reliable_alleles_text$COUNT2 >= 0, 0, 1)
      p2 <- p2 + geom_text(data = non_reliable_alleles_text, aes_string(label = "text", x = "GENE", hjust = "hjust"), angle = 0, size = 2.25)
      short_reads = T
    }
    ########################################################################################################
    
    ### Plot All panels
    
    if (html_output) {
      options(warn = -1)
      ## Prepare panels for html plot
      
      p = p + theme(axis.title.x = element_blank())
      p.l <- ggplotly(p, height = 800, width = 400, tooltip = "text") %>% plotly::layout(showlegend = FALSE)
      
      pk = pk + theme(axis.title = element_blank())
      pk <- pk + scale_fill_manual(name = "log<sub>10</sub>(lK)", values = c('#FFFFFF',RColorBrewer::brewer.pal(9,'Blues')), drop = FALSE)
      pk.l <- ggplotly(pk, height = 800, width = 400, tooltip = "text") %>% plotly::layout(showlegend = TRUE)
      pk.l$x$layout$annotations[[1]]$text = p.l$x$layout$annotations[[1]]$text
      pk.l$x$layout$annotations[[2]]$text = p.l$x$layout$annotations[[2]]$text
      
      p2 <- p2 + ylab("log<sub>10</sub>(Count+1)")
      p2.l <- ggplotly(p2, height = 1000, width = 700, tooltip = "text") %>% plotly::layout(margin = list(b = 50), yaxis = list(title = paste0(c(rep("&nbsp;", 3), "Gene",
                                                                                                                                                 rep("&nbsp;", 3), rep("\n&nbsp;", 1)), collapse = "")), showlegend = TRUE)
      
      p2.l$x$layout$xaxis$ticktext = c(lapply(p2.l$x$layout$xaxis$ticktext[1:match("0", p2.l$x$layout$xaxis$ticktext) - 1], function(x) paste0("-",
                                                                                                                                               x)), p2.l$x$layout$xaxis$ticktext[match("0", p2.l$x$layout$xaxis$ticktext):length(p2.l$x$layout$xaxis$ticktext)])
      
      p2.l$x$data[[grep("dot",p2.l$x$data)]]$y[2] = p2.l$x$layout$yaxis$range[2]
      
      add_border <- function(p_l){
        p_l$x$attrs <- lapply(p_l$x$attrs,
                              function(x){
                                x$mode <- "markers"
                                x
                              })
        
        for (i in 1:length(p_l$x$data)) {
          p_l$x$data[[i]]$marker$line$width <- suppressMessages(0.2)
          p_l$x$data[[i]]$marker$line$color <- ifelse(p_l$x$data[[i]]$marker$line$color=='transparent', 'black', p_l$x$data[[i]]$marker$line$color)
        }
        return(p_l)
      }
      
      p.l <- add_border(p.l)
      pk.l <- add_border(pk.l)
      mgsub <- function(pattern, replacement, x, ...) {
        if (length(pattern) != length(replacement)) {
          stop("pattern and replacement do not have the same length.")
        }
        result <- x
        for (i in 1:length(pattern)) {
          result <- gsub(pattern[i], replacement[i], result, ...)
        }
        result
      }
      
      
      text_for_hovertext <- function(labels, count.df) {
        for (i in 1:length(labels)) {
          label <- labels[i]
          gene <- strsplit(strsplit(label, "<")[[1]][1], " ")[[1]][2]
          allele <- strsplit(label, "Allele: ")[[1]][2]
          if (!is.na(NA)) {
            count <- strsplit(strsplit(label, "<br />Count: ")[[1]][2], "<")[[1]][1]
            if (count == "NA")
              next
            count <- as.numeric(count)
            if (count%%1 != 0)
              count <- count.df %>% filter(.data$GENE == gene & .data$ALLELES == allele & round(.data$COUNT3, nchar(as.character(count)) - 2) == count) %>% select(.data$COUNT) else count <- count.df %>% filter(.data$GENE == gene & .data$ALLELES == allele & .data$COUNT3 == count) %>% select(.data$COUNT)
              labels[i] <- paste0("Gene: ", gene, "<br />Allele: ", allele, "<br />Count: ", count[1, ])
          }
        }
        return(labels)
      }
      
      text_for_hovertext_non_reliable <- function(labels, text, annot) {
        for (i in 1:length(labels)) {
          label <- labels[i]
          label.tmp <- toupper(label)
          gene <- strsplit(strsplit(label.tmp, "GENE: ")[[1]][2], '[<]')[[1]][1]
          allele <- gsub("text: ","",strsplit(grep(gene,text,value=T), "<")[[1]][1])
          
          labels[i] <- paste0("Gene: ", gene, "<br />Allele: ", annot[allele])
        }
        return(labels)
      }
      
      count.df$COUNT3 <- abs(count.df$COUNT2)
      
      change_hovertext <- function(pp,annot){
        ind <- grep('text: ',pp$x$data)
        for(i in ind){
          pp$x$data[[i]]$type <-  suppressWarnings("text")
          pp$x$data[[i]]$hoveron <- "fill"
          pp$x$data[[i]]$hoverinfo <- "skip"
        }
        return(pp)
      }
      
      bottom_annot <- c()
      if(short_reads){
        bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
        names(bottom_annot) <- unique(non_reliable_alleles_text$text)
        p.l <- change_hovertext(p.l,bottom_annot)
        p2.l <- change_hovertext(p2.l,bottom_annot)
      }
      
      # for (i in 1:length(p2.l$x$data)) {
      #   p2.l$x$data[[i]]$text <- mgsub(c("border: ","white|black","^[<]br [/][>]","GENE", "COUNT2", "ALLELES","factor[(]ALLELES[,] levels [=] AlleleCol[)]", "yintercept: 0"),
      #                                  c("","","","Gene", "Count","Allele","Allele", ""), p2.l$x$data[[i]]$text)
      #   p2.l$x$data[[i]]$text <- text_for_hovertext(p2.l$x$data[[i]]$text, count.df)
      # }
      #
      # for (i in 1:length(p.l$x$data)) {
      #     p.l$x$data[[i]]$text <- mgsub(c("count\\: [0-9]\\.[0-9]","^[<]br [/][>]","GENE", "ALLELES","factor[(]ALLELES[,] levels [=] AlleleCol[)]"),
      #                                   c("","","Gene", "Allele", "Allele"), p.l$x$data[[i]]$text)
      # }
      #
      # for (i in 1:length(pk.l$x$data)) {
      #     pk.l$x$data[[i]]$text <- mgsub(c("count\\: [0-9]","^[<]br [/][>]","GENE", "K[_]GROUPED"),
      #                                      c("","","Gene", "K"), pk.l$x$data[[i]]$text)
      # }
      
      
      p.l.c <- suppressWarnings(subplot(p2.l, p.l, pk.l, widths = c(0.4, 0.2, 0.2), shareY = T, titleX = TRUE, margin = 0.01, which_layout = 1) %>%
                                  plotly::layout(title = sample_name, titlefont=list(size=16), margin = list(t = 45)))
      
      
      for (i in 1:length(p2.l$x$data)) {
        p.l.c$x$data[[i]]$showlegend <- FALSE
      }
      
      for(i in 1:(length(p2.l$x$data)+length(p.l$x$data))){
        p.l.c$x$data[[i]]$name <- gsub("[()|,]|white|black","",p.l.c$x$data[[i]]$name)
      }
      
      p.l.c$x$layout$annotations[[6]]$text = "log<sub>10</sub>(lK)----"
      p.l.c$x$layout$annotations[[6]]$xanchor = "center"
      p.l.c$x$layout$annotations[[6]]$legendtitle = TRUE
      p.l.c$x$layout$annotations[[6]]$y = 1 - 0.033 * (length(grep('[[]',grep('TRUE',p.l.c$x$data,value = T),invert = T))) #(length(AlleleCol) + 2.4)  #0.52
      p.l.c$x$layout$annotations[[6]]$x = 0.985
      p.l.c$x$layout$annotations[[6]]$xref = 'paper'
      p.l.c$x$layout$annotations[[6]]$yref = 'paper'
      p.l.c$x$layout$annotations[[6]]$font$size = 16
      
      
      p.l.c$x$layout$annotations[[3]] <- p.l.c$x$layout$annotations[[6]]
      p.l.c$x$layout$annotations[[3]]$text = "Alleles----"
      p.l.c$x$layout$annotations[[3]]$y = 0.99
      p.l.c$x$layout$annotations[[3]]$x = 0.98
      p.l.c$x$layout$annotations[[3]]$legendTitle = FALSE
      p.l.c$x$layout$annotations[[3]]$font$size = 16
      
      if(short_reads){
        bottom_annot_collapsed <- ifelse(length(bottom_annot) > 1,
                                         paste0(sapply(split(bottom_annot, ceiling(seq_along(bottom_annot)/1)),
                                                       function(x) paste0(x,collapse = '\t')),collapse = '\n'),
                                         paste0(bottom_annot,collapse = '\t'))
        
        p.l.c$x$layout$annotations[[7]] <- p.l.c$x$layout$annotations[[6]]
        p.l.c$x$layout$annotations[[7]]$text = bottom_annot_collapsed
        p.l.c$x$layout$annotations[[7]]$y = p.l.c$x$layout$annotations[[6]]$y - 0.05 * (length(grep('[[]',grep('TRUE',p.l.c$x$data,value = T),invert = T)))
        p.l.c$x$layout$annotations[[7]]$x = 1.1
        p.l.c$x$layout$annotations[[7]]$legendTitle = FALSE
        p.l.c$x$layout$annotations[[7]]$font$size = 12
      }
      
      
      plot_list[[sample_name]] <- p.l.c
      
    } else {
      p.legend <- get_legend(p + theme(legend.key = element_rect("#DCDCDC")))
      p = p + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.y = element_blank())
      
      pk = pk + scale_fill_manual(name = expression("log"[10] * "(lK)"), values = c('#FFFFFF',RColorBrewer::brewer.pal(9,'Blues')), drop = FALSE)
      pk.legend <- get_legend(pk + theme(legend.key = element_rect("gray")))
      pk = pk + theme(legend.position = "none",axis.title = element_blank())
      
      p.legends <- plot_grid(pk.legend, p.legend, ncol = 1, align = "hv")
      
      p1 <- plot_grid(p2, p, pk, nrow = 1, rel_widths = c(0.2, 0.125, 0.1), align = "hv", axis = "b")
      
      p <- plot_grid(p1, p.legends, ncol = 2, rel_widths = c(0.9, 0.1))
      
      if(short_reads){
        bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
        # Create text for annotating the short reads labels
        lab <- grid::textGrob(ifelse(length(bottom_annot) > 10,
                                     paste0(sapply(split(bottom_annot, ceiling(seq_along(bottom_annot)/10)),
                                                   function(x) paste0(x,collapse = '\t')),collapse = '\n'),
                                     paste0(bottom_annot,collapse = '\t')),
                              x = unit(.1, "npc"), just = c("left"),
                              gp = grid::gpar(fontsize = 10, col = "black"))
        
        gp <- ggplotGrob(p)
        # Add a row below the 2nd from the bottom
        gp <- gtable::gtable_add_rows(gp, unit(2, "grobheight", lab), -2)
        
        # Add 'lab' grob to that row, under the plot panel
        gp <- gtable::gtable_add_grob(gp, lab, t = -2, l = gp$layout[gp$layout$name == "panel",]$l)
        
        plot_list[[sample_name]] <- gp
      } else{
        plot_list[[sample_name]] <- p
      }
      
    }
    
  }
  if (length(plot_list) != 1) {
    
    
    if(!missing(dir)){
      dir <- file.path(dir, "haplotype_output")
      dir.create(dir)
      
    }else{
      dir <- tempdir()
    }
    
    if (html_output) {
      
      for (sample_name in names(plot_list)) {
        htmlwidgets::saveWidget(plot_list[[sample_name]], paste0(dir, '/', sample_name, ".html"), selfcontained = F)
      }
    } else {
      pdf(paste0(dir, "/haplotype_output.pdf"), height = 20, width = 15)
      for (sample_name in names(plot_list)) {
        
        title <- ggdraw() + draw_label(sample_name, fontface = "bold")
        plot(plot_grid(title, plot_list[[sample_name]], ncol = 1, rel_heights = c(0.05, 1)))
      }
      
      dev.off()
    }
    
  } else if (html_output){
    return(plot_list[[1]])
  }else{
    title <- ggdraw() + draw_label(sample_name, fontface = "bold")
    plot(plot_grid(title, plot_list[[1]], ncol = 1, rel_heights = c(0.05, 1)))
  }
}

########################################################################################################
#' Graphical output of alleles division by chromosome
#'
#' The \code{hapHeatmap} function generates a graphical output of the alleles per gene in multiple samples.
#'
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    lk_cutoff            the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.
#' @param    mark_low_lk          if TRUE, a texture is add for low lK values. Defualt is TRUE.
#' @param    size_annot           size of bottom annotation text. Defualt is 1.5 .
#' @param    color_y              named list of the colors for y axis labels.
#' @param    order_subject        order subject by a vecor.
#' @param    file                 file path for rendering the plot to pdf. If non is supplied than the plot is retured as object. Defualt is NULL.
#'
#' @return
#'
#' A list with the following:
#'
#' \itemize{
#'   \item \code{'p'}:        heat-map visualization of the haplotype inference for multiple samples.
#'   \item \code{'width'}:    Optimal width value for rendering plot.
#'   \item \code{'height'}:   Optimal width value for rendering plot.
#' }
#'
#' When a file is supplied the graph is also rendered to pdf.
#'
#' @details
#'
#' A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @examples
#' # Plotting haplotpe heatmap
#' p <- hapHeatmap(samplesHaplotype)
#'
#' cowplot::ggdraw(p$p)
#' @export
hapHeatmap <- function(hap_table, chain = c("IGH", "IGK", "IGL", "TRB"), gene_sort = "position", removeIGH = TRUE, lk_cutoff = 1, mark_low_lk = TRUE, size_annot = 1.5, color_y = NULL, order_subject = NULL , file = NULL) {
  
  
  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)
  
  lk_cutoff = as.numeric(lk_cutoff)
  
  id_GENE_col <- which(names(hap_table)=="GENE")
  hapBy_col_id <- c(id_GENE_col+1,id_GENE_col+2)
  hapBy_cols <- names(hap_table)[hapBy_col_id]
  hapBy_alleles <- gsub("_", "*", hapBy_cols)
  samples <- unique(hap_table$SUBJECT)
  
  k_assign <- function(x){
    # getting the lk value for each allele
    # input haplotype row
    # output lk value
    if(x[5] %fin% c("Unk","Del","NR")){
      k_m <- paste0('K',1:4)
    }else{
      k_m <- paste0('K',which(x[-c(1:10)]==x[5],arr.ind = T))
    }
    k <- min(as.numeric(x[k_m]), na.rm = T)
  }
  
  # sort the data
  panels <- sortDFByGene(hap_table, chain = chain, method = gene_sort, removeIGH = removeIGH, geno = T, peseudo_remove = T)
  panels$GENE <- factor(panels$GENE, levels = gsub(chain, "", GENE.loc[[chain]]))
  # rename genes to numbers
  gene_loc <- 1:length(unique(panels$GENE)[order(match(unique(panels$GENE), levels(panels$GENE)))])
  names(gene_loc) <- unique(panels$GENE)[order(match(unique(panels$GENE), levels(panels$GENE)))]
  panels$GENE_LOC <- gene_loc[as.character(panels$GENE)]
  # fix na in K columns
  panels[is.na(panels)] <- Inf
  # melt haplotype columns to one
  panels <- melt(panels, measure.vars = hapBy_cols, variable.name = "hapBy", value.name = 'hapBy_alleles')
  
  # separating the allele column
  panels <- splitstackshape::cSplit(panels[,c('SUBJECT','GENE','GENE_LOC',"hapBy",'hapBy_alleles',"ALLELES",paste0("K",1:4))],
                                    'ALLELES', sep = ",", fixed = T, type.convert = F, drop = F)
  # separating the hap alleles
  panels <- splitstackshape::cSplit(panels, 'hapBy_alleles', sep = ",", direction = "long", fixed = T, type.convert = F)
  # getting the lk value
  panels[,"K":=apply(panels,1,k_assign)]
  # clean lk inf
  invisible(lapply(names(panels),function(.name) data.table::set(panels, which(is.infinite(panels[[.name]])), j = .name,value=NA)))
  # remove unnecessary columns
  cols <- c("SUBJECT","GENE","GENE_LOC","hapBy","hapBy_alleles","K")
  panels <- panels[, .SD, .SDcols= cols]

  ######sort the heatmap for plotting
  panels_m <- panels[, "n":=  .N, by = c("SUBJECT", "GENE", "hapBy")][] # count number of alleles for group
  panels_m$ALLELES_G <- panels_m$hapBy_alleles # for grouping
  panels_m$text <- ''
  panels_m$text_bottom <- panels_m$hapBy_alleles
  # change ambiguous hapBy_alleles call
  id_nra <- grepl("[0-9][0-9]_[0-9][0-9]$", panels_m$hapBy_alleles)
  nra <- F
  if (any(id_nra)) {
    # number ambiguous hapBy_alleles
    num_text <- paste0('[*',1:length(unique(panels_m$hapBy_alleles[id_nra])),']')
    names(num_text) <- unique(panels_m$hapBy_alleles[id_nra])
    # text for plot
    panels_m$text[id_nra] <- num_text[panels_m$hapBy_alleles[id_nra]]
    # text for legend
    panels_m$text_bottom[id_nra] <- paste(num_text[panels_m$hapBy_alleles[id_nra]],panels_m$hapBy_alleles[id_nra])
    # change allele to NRA - non reliable allele
    panels_m$hapBy_alleles[id_nra] <- "NRA"
    # indicates that nra exists
    nra <- T
  }
  # create allele palette
  allele_palette <- alleleHapPalette(panels_m$hapBy_alleles)
  
  # sort novel allele calls for plot
  val_novel <- grep('^[0-9]+[_][A-Z]|^[0-9]+[_][0-9]+[A-Z]',panels_m$hapBy_alleles, value = T)
  novel <- F
  novel_allele_text <- c()
  novel_symbol <- "\u005E"
  if(length(val_novel)!=0){
    # sort the palettle colors for novel alleles
    id <- grep('^[0-9]+[_][A-Z]|^[0-9]+[_][0-9]+[A-Z]',names(allele_palette$transper))
    allele_palette$transper[id] <- 1
    # cerate code index for novel allele
    code_allele <- paste0(novel_symbol,1:length(id))
    names(code_allele) <-allele_palette$AlleleCol[id]
    new_allele <- paste0(novel_symbol,1:length(id),'-',allele_palette$AlleleCol[id])
    names(new_allele) <-allele_palette$AlleleCol[id]
    # change the text for plot
    ids <- panels_m$hapBy_alleles %fin% names(new_allele)
    rep <- new_allele[panels_m$hapBy_alleles[ids]]
    rep2 <- code_allele[panels_m$hapBy_alleles[ids]]
    # add new allele code to data
    panels_m[ids, c("hapBy_alleles","text_bottom","text") := list(rep,rep,rep2)]
    # change annotation in legend colors
    allele_palette$AlleleCol[id] <- new_allele
    names(allele_palette$transper)[id] <- new_allele
    # indicates that novel exists
    novel <- T
  }
  
  panels_m$hapBy_alleles <- factor(panels_m$hapBy_alleles, levels = allele_palette$AlleleCol)
  
  # samples names and number
  samples <- unique(panels_m$SUBJECT)
  samples_n <- length(samples)
  # genes names and number
  genes <- unique(panels_m$GENE)
  genes_n <- length(genes)

  # order the data by gene loc
  setorderv(panels_m, c("SUBJECT","GENE_LOC"))
  setkey(panels_m, "SUBJECT")
  
  # sort data for matrix
  panels_m[,"line":=12/panels_m$n]
  allele_code <- 1:length(allele_palette$AlleleCol)
  names(allele_code) <- gsub(paste0("\\^","[0-9]+[-]"),"",allele_palette$AlleleCol)
  # sort the alleles in gene box
  panels_m[,"A_CODE":=allele_code[hapBy_alleles]+1]
  panels_m[grep("[0-9][0-9]_[0-9][0-9]$",panels_m$hapBy_alleles,perl = T),"A_CODE":=allele_code["NRA"]]
  setorderv(panels_m, c("SUBJECT","GENE_LOC","A_CODE"))
  
  # duplicate the data by 12 box to gene
  panels_m[,"id" := 1:.N, by = c("SUBJECT", "GENE", "hapBy")]
  panels_f = panels_m[,c("n_line" = 1:get("line")), by = c("SUBJECT", "hapBy", "GENE", "GENE_LOC", "ALLELES_G", "A_CODE", "text_bottom", "K"), nomatch = 0]
  
  if(!is.null(order_subject)) panels_f <- panels_f[order(match(panels_f$SUBJECT, order_subject))]
  
  
  # transform allele codes to matrix, 12 box for each gene. each row is an individual
  upper_m <- matrix(panels_f[panels_f$hapBy==hapBy_cols[1]][[6]],ncol = 12*genes_n,byrow = T,
                    dimnames = list(unique(panels_f[panels_f$hapBy==hapBy_cols[1]][[1]]),panels_f[panels_f$hapBy==hapBy_cols[1]][[3]][1:(12*genes_n)]))
  lower_m <- matrix(panels_f[panels_f$hapBy==hapBy_cols[2]][[6]],ncol = 12*genes_n,byrow = T,
                    dimnames = list(unique(panels_f[panels_f$hapBy==hapBy_cols[2]][[1]]),panels_f[panels_f$hapBy==hapBy_cols[2]][[3]][1:(12*genes_n)]))
  
  #### sort legend
  # fit width of legend column text by the longest allele
  longest_allele <- max(nchar(allele_palette$AlleleCol))*3+40
  # legend length
  leg_length <- next_divisor(longest_allele,genes_n*12)
  # create legend values for matrix
  leg <- function(x,allele_text_length) c(rep(x,4),rep(0,allele_text_length))
  seqs <- lapply(allele_code+1,leg, allele_text_length = leg_length-4)
  # add values for legend to complete the matrix, white boxes
  add = ceiling(length(unlist(seqs))/(genes_n*12))*(genes_n*12) - length(unlist(seqs))
  add = ifelse(add<0,0,add)
  # legend matrix
  m2 <- matrix(c(unlist(seqs),rep(0,add)), ncol = genes_n*12, byrow = T)
  
  # start and end values for plot parts
  start <- c(0.01)
  end <- c(0.15,0.13)
  size = 1
  short_reads_rows = 0
  ## add short read text annotation at the bottom
  if(nra){
    bottom_annot <- unique(grep("[0-9][0-9]_[0-9][0-9]$", panels_f$text_bottom , value = T, perl = T))
    
    # Create text for annotating the short reads labels
    annot <- splitlines(bottom_annot, genes_n*12/(4*size_annot))
    annot <- annot[!is.na(dplyr::na_if(annot,""))]
    x_annot = max(nchar(annot))/(12*genes_n)
    short_reads_rows = length(annot)
    annot <- paste0(annot,collapse = '\n')
    # add the start and value for the third part
    end <- c(0.2, 0.22, 0.06)
    start <- c(0.03,0.01)
  }
  
  # set the height and width of plot
  height <- samples_n * 0.1 + 10 + nrow(m2)*0.2 + short_reads_rows*0.4 # number of samples, number of rows in legend, number of rows in bottom annotation
  width <- genes_n * 0.3 + 5 # numer of genes
  size_text = nrow(upper_m)/(height*width)+0.5 # text size for heatmap annoations
  size_text_leg = ncol(m2)/(width*longest_allele)+1 # text size for legend annotations
  
  if(!is.null(file)) pdf(file, onefile = F, width = width, height = height, family = "serif")
  # plot layout
  layout.matrix <- matrix(c(1,2, 3, 4), nrow = 4, ncol = 1)
  graphics::layout(mat = layout.matrix,
                   heights = c(2, 2, 1, 1) # Heights of the three rows
  )
  # heatmap upper plot
  par(mar=c(2,6,8,6))
  image(t(upper_m),col = names(allele_palette$AlleleCol), breaks = 1:(length(allele_palette$AlleleCol)+1),axes=F)
  # add grid lines for genes
  grid(lwd=1,nx = genes_n,ny=0,col = "white",lty = 1)
  # add axis annotations
  axis(3,(0:(genes_n-1))/genes_n+6/(12*genes_n),names(gene_loc),las=3) # top
  axis(1,(0:(genes_n-1))/genes_n+6/(12*genes_n),names(gene_loc),las=3) # bottom
  title(gsub('_','*',hapBy_cols[1]), adj = 0.5)
  # color y tick labels if supplied
  colors <- "black"
  if(!is.null(color_y)) colors <- color_y[rownames(upper_m)]
  Map(axis, side=2, at=(0:(samples_n-1))/(samples_n-1), col.axis=colors, labels=rownames(upper_m), lwd=0, las=1, cex.axis=0.8) #left
  axis(2,at=(0:(samples_n-1))/(samples_n-1),labels=FALSE)
  
  
  # draw lines for low lk values
  sub_geno = panels_m[panels_m$hapBy==hapBy_cols[1] & panels_m$K<lk_cutoff,]
  NR = samples_n
  NC = genes_n*12
  apply(sub_geno, 1,function(x){
    
    I = which(x["SUBJECT"]==samples)-1    # row index
    J = (as.numeric(x["GENE_LOC"])-1)*12            # column index
    draw_segment(NR,NC,I,J,lwd=1,col="white")}
  )
  
  # ad text annotations
  ids_text <- !grepl('^[0-9]|Del|Unk',panels_m$text_bottom)
  sub_geno = panels_m[panels_m$hapBy==hapBy_cols[1] & ids_text,]
  
  NR = samples_n
  NC = genes_n*12
  apply(sub_geno, 1,function(x){
    
    I = which(x["SUBJECT"]==samples)-1    # row index
    J = (as.numeric(x["GENE_LOC"])-1)*12            # column index
    ALLELE =  as.numeric(x["id"])                   # allele index
    N_ALLELES = as.numeric(x["n"])                  # number of alleles
    TEXT =  x["text"]                   # text
    Write_text(NR,NC,I,J,ALLELE,N_ALLELES,TEXT,cex=size_text)}
  )
  
  # heatmap lower plot
  par(mar=c(2,6,8,6))
  image(t(lower_m),col = names(allele_palette$AlleleCol), breaks = 1:(length(allele_palette$AlleleCol)+1),axes=F)
  # add grid lines for genes
  grid(lwd=1,nx = genes_n,ny=0,col = "white",lty = 1)
  # add axis annotations
  axis(3,(0:(genes_n-1))/genes_n+6/(12*genes_n),names(gene_loc),las=3) # top
  axis(1,(0:(genes_n-1))/genes_n+6/(12*genes_n),names(gene_loc),las=3) # bottom
  title(gsub('_','*',hapBy_cols[2]), adj = 0.5)
  # color y tick labels if supplied
  colors <- "black"
  if(!is.null(color_y)) colors <- color_y[rownames(lower_m)]
  Map(axis, side=2, at=(0:(samples_n-1))/(samples_n-1), col.axis=colors, labels=rownames(lower_m), lwd=0, las=1, cex.axis=0.8) #left
  axis(2,at=(0:(samples_n-1))/(samples_n-1),labels=FALSE)
  
  # draw lines for low lk values
  sub_geno = panels_m[panels_m$hapBy==hapBy_cols[2] & panels_m$K<lk_cutoff,]
  NR = samples_n
  NC = genes_n*12
  apply(sub_geno, 1,function(x){
    
    I = which(x["SUBJECT"]==samples)-1    # row index
    J = (as.numeric(x["GENE_LOC"])-1)*12            # column index
    draw_segment(NR,NC,I,J,lwd=1,col="white")}
  )
  
  # ad text annotations
  sub_geno = panels_m[panels_m$hapBy==hapBy_cols[2] & ids_text,]
  NR = samples_n
  NC = genes_n*12
  apply(sub_geno, 1,function(x){
    
    I = which(x["SUBJECT"]==samples)-1    # row index
    J = (as.numeric(x["GENE_LOC"])-1)*12            # column index
    ALLELE =  as.numeric(x["id"])                   # allele index
    N_ALLELES = as.numeric(x["n"])                  # number of alleles
    TEXT =  x["text"]                   # text
    Write_text(NR,NC,I,J,ALLELE,N_ALLELES,TEXT,cex=size_text)}
  )
  
  # legend plot
  par(mar=c(3,6,2,6))
  image(t(m2),col = names(allele_palette$AlleleCol), breaks = 1:(length(allele_palette$AlleleCol)+1),axes=F)
  # add grid lines for each row
  grid(lwd=1,ny = nrow(m2), nx = 0,col = "black",lty = 2)
  # add text anotation for legend
  NR = nrow(m2)
  NC = genes_n*12
  names(allele_code) <- allele_palette$AlleleCol
  invisible(tapply(allele_code,names(allele_code),function(x){
    ii = which(m2==x+1,arr.ind = T)[1,]
    I = ii[[1]]-1              # row index
    J = ii[[2]]                # column index
    ALLELE =  1                # allele index
    N_ALLELES = 1              # number of alleles
    TEXT =  names(x)            # text
    STEP_X<-1/(NC-1)
    STEP_Y<-ifelse(1/(NR-1)==Inf,0,1/(NR-1))
    text(STEP_X*J+STEP_X*leg_length*0.5,
         STEP_Y*I,
         TEXT,cex = size_text_leg)
  }
  ))
  # bottom text for nra, only if exists
  if(nra){
    # add the text at the bottom
    par(mar=c(1,6,3,6))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = x_annot, y = 0.5, annot, cex = size_annot, col = "black", pos = 4)
  }
  p1.base <- recordPlot()
  invisible(dev.off())
  return(list(p = p1.base, width = width, height = height))
  if(!is.null(file)) dev.off()
  # embed the fonts to file
  if(!is.null(file)) embedFonts(file)
}
