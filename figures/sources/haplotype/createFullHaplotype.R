
##########################################################################
#' Anchor gene haplotype inference
#'
#' The \code{createFullHaplotype} functions infers haplotype based on an anchor gene.
#'
#'
#' @param    clip_db               a \code{data.frame} in Change-O format. See details.
#' @param    toHap_col             a vector of column names for which a haplotype should be inferred. Default is V_CALL and D_CALL.
#' @param    hapBy_col             column name of the anchor gene. Default is J_CALL.
#' @param    hapBy                 a string of the anchor gene name. Default is IGHJ6.
#' @param    toHap_GERM            a vector of named nucleotide germline sequences matching the allele calls in \code{toHap_col} columns in clip_db.
#' @param    relative_freq_priors  if TRUE, the priors for Bayesian inference are estimated from the relative frequencies in clip_db. Else, priors are set to \code{c(0.5,0.5)}. Defualt is TRUE
#' @param    kThreshDel            the minimum lK (log10 of the Bayes factor) to call a deletion. Defualt is 3.
#' @param    rmPseudo              if TRUE non-functional and pseudo genes are removed. Defualt is TRUE.
#' @param    deleted_genes         double chromosome deletion summary table. A \code{data.frame} created by \code{deletionsByBinom}.
#' @param    nonReliable_Vgenes    a list of known non reliable gene assignments. A \code{list} created by \code{nonReliableVGenes}.
#' @param    min_minor_fraction    the minimum minor allele fraction to be used as an anchor gene. Default is 0.3
#' @param    chain                 the IG chain: IGH,IGK,IGL. Default is IGH.
#'
#' @return
#' A \code{data.frame}, in which each row is the haplotype inference summary of a gene from the column selected in \code{toHap_col}.
#'
#'The output containes the following columns:
#' \itemize{
#'  \item \code{SUBJECT}:        the subject name.
#'  \item \code{GENE}:           the gene name.
#'  \item Anchor gene allele 1:  the haplotype inference for chromosome one. The column name is the anchor gene with the first allele.
#'  \item Anchor gene allele 2:  the haplotype inference for chromosome two. The column name is the anchor gene with the second allele.
#'  \item \code{ALLELES}:        allele calls for the gene.
#'  \item \code{PRIORS_ROW}:     priors based on relative allele usage of the anchor gene.
#'  \item \code{PRIORS_COL}:     priors based on relative allele usage of the inferred gene.
#'  \item \code{COUNTS1}:        the appereance count on each chromosome of the first allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K1}:             the Bayesian factor value for the first allele (from \code{ALLELES}) inference.
#'  \item \code{COUNTS2}:        the appereance count on each chromosome of the second allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K2}:             the Bayesian factor value for the second allele (from \code{ALLELES}) inference.
#'  \item \code{COUNTS3}:        the appereance count on each chromosome of the third allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K3}:             the Bayesian factor value for the third allele (from \code{ALLELES}) inference.
#'  \item \code{COUNTS4}:        the appereance count on each chromosome of the fourth allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K4}:             the Bayesian factor value for the fourth allele (from \code{ALLELES}) inference.
#'}
#'
#' @details
#' Function accepts a \code{data.frame} in Change-O format (\url{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{'SUBJECT'}: The subject name
#'   \item \code{'V_CALL'}: V allele call(s) (in an IMGT format)
#'   \item \code{'D_CALL'}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{'J_CALL'}: J allele call(s) (in an IMGT format)
#' }
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db, HVGERM, HDGERM)
#'
#' # Selecting a single individual
#' clip_db = samples_db[samples_db$SUBJECT=='I5', ]
#'
#' # Infering haplotype
#' haplo_db = createFullHaplotype(clip_db,toHap_col=c('V_CALL','D_CALL'),
#' hapBy_col='J_CALL',hapBy='IGHJ6',toHap_GERM=c(HVGERM,HDGERM))
#'
#'
#' @export

createFullHaplotype <- function(clip_db, toHap_col = c("V_CALL", "D_CALL"), hapBy_col = "J_CALL", hapBy = "IGHJ6", toHap_GERM, relative_freq_priors = TRUE,
                                kThreshDel = 3, rmPseudo = TRUE, deleted_genes = c(), nonReliable_Vgenes = c(), min_minor_fraction = 0.3, chain = c("IGH", "IGK", "IGL", "TRB")) {
  
  
  # Check if germline was inputed
  if (missing(toHap_GERM))
    stop("Missing toHap_GERM, please input germline sequences")
  
  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)
  
  if (!("SUBJECT" %in% names(clip_db))) {
    
    clip_db$SUBJECT <- "S1"
  }
  
  haplo_db <- c()
  clip_db <- clip_db %>% select(.data$SUBJECT, !!eval(c(hapBy_col,toHap_col)))
  for (sample_name in unique(clip_db$SUBJECT)) {
    
    if (is.list(nonReliable_Vgenes)){
      nonReliable_Vgenes_vec <- nonReliable_Vgenes[[sample_name]]
    } else nonReliable_Vgenes_vec <- nonReliable_Vgenes
    
    if (is.data.frame(deleted_genes)) {
      deleted_genes_vec <- deleted_genes %>% filter(.data$SUBJECT == sample_name, .data$DELETION == "Deletion") %>% select(.data$GENE) %>% pull()
      if (is.null(nonReliable_Vgenes_vec))
        nonReliable_Vgenes_vec <- deleted_genes %>% filter(.data$SUBJECT == sample_name, .data$DELETION == "Non reliable") %>% select(.data$GENE) %>% pull()
    } else deleted_genes_vec <- c()
    
    
    ### Check if haplotype can be infered by the specific gene in the data set.  Only relevant genes with one assignment
    clip_db_sub <- clip_db %>% filter(.data$SUBJECT==sample_name, !grepl(',', !!as.name(hapBy_col), perl = T))
    hapBy_priors <- clip_db_sub %>% filter(grepl(paste0(hapBy, "\\*"), !!as.name(hapBy_col), perl = T)) %>% dplyr::count(!!as.name(hapBy_col)) %>% mutate(freq = n / sum(n)) %>% select(-n)
    hapBy_priors <- setNames(hapBy_priors[[2]], hapBy_priors[[1]])
    print(hapBy_priors)
    hapBy_alleles <- names(hapBy_priors)
    if (length(hapBy_alleles) != 2)
      stop("Can not haplotype by more or less than two alleles")
    if (min(hapBy_priors) < min_minor_fraction)
      stop("Can not haplotype, minor allele fraction lower than the cutoff set by the user")
    
    GENES <- unique(gsub("\\*.*", "\\1", grep("^(?=.*IG)|^(?=.*TR)", unique(unlist(clip_db_sub[clip_db_sub[, hapBy_col] %in% hapBy_alleles,toHap_col], use.names = F)), value = T, perl = T), perl = T))
    
    GENES.ref <- unique(gsub("\\*.*", "\\1", names(toHap_GERM), perl = T))
    if (rmPseudo) {
      GENES <- GENES[!grepl("OR|NL", GENES)]
      GENES <- GENES[!(GENES %in% PSEUDO[[chain]])]
      
      GENES.ref <- GENES.ref[!grepl("OR|NL", GENES.ref)]
      GENES.ref <- GENES.ref[!(GENES.ref %in% PSEUDO[[chain]])]
    }
    GENES.df.num <- data.table::rbindlist(lapply(intersect(GENES, GENES.ref),function(G){
      
      if(G %in% deleted_genes_vec || G %in% nonReliable_Vgenes_vec){
        
        relFreqDf.tmp <- data.frame(matrix(c(sample_name, G,
                                             rep(ifelse(G %in% nonReliable_Vgenes_vec, 'NR',
                                                        ifelse(G %in% deleted_genes_vec, 'Del')), 2),
                                             rep(NA, 11)),nrow = 1),stringsAsFactors = F)
        
        relFreqDf.tmp <- asNum(relFreqDf.tmp)
        return(relFreqDf.tmp)
      } else{
        
        toHap_col_tmp <-   toHap_col[stringi::stri_detect_fixed(pattern = substr(G,4,4),str = toHap_col)]
        
        clip_db_sub.G <- clip_db_sub %>% filter(grepl(paste0("^(",G,"\\*[[:digit:]]*[\\_[[:alnum:]]*]*,?)+$"),!!as.name(toHap_col_tmp),perl = T))
        if(substr(G,4,4)=='V'){
          tmp <- clip_db_sub.G %>% filter(stringi::stri_detect_regex(pattern = ",",str = !!as.name(toHap_col_tmp), negate = T))
          tmp2 <- data.table(clip_db_sub.G %>% filter(stringi::stri_detect_regex(pattern = ",",str = !!as.name(toHap_col_tmp))), key=toHap_col_tmp)
          tmp2 <- tmp2[, n:=.N, by=eval(toHap_col_tmp)][,"prop" := n/nrow(clip_db_sub.G)][get("prop") > 0.4,]
          if(length(tmp2[[eval(toHap_col_tmp)]])>0) tmp2[[toHap_col_tmp]] <- alleleCollapse(tmp2[[eval(toHap_col_tmp)]])
          clip_db_sub.G <- rbind(tmp,as.data.frame(tmp2[, c("n", "prop"):=NULL]))
          
        } else clip_db_sub.G <-  clip_db_sub.G %>% filter(!grepl(',', !!as.name(toHap_col_tmp)))
        
        tmp <- clip_db_sub.G %>% filter(grepl(paste0('^(?=.*',hapBy,'\\*)'), !!as.name(hapBy_col), perl = T)) %>% select(!!as.name(toHap_col_tmp), !!as.name(hapBy_col)) %>% table()
        
        if(nrow(tmp)==0){return()
        }else{
          # if one column add the second
          if (ncol(tmp) == 1) {
            toadd <- setdiff(hapBy_alleles, colnames(tmp))
            tmp <- cbind(tmp, rep(0, nrow(tmp)))
            colnames(tmp)[2] <- toadd
            
            tmp <- as.data.frame(tmp)
            tmp <- tmp[order(colnames(tmp))]
            tmp <- as.matrix(tmp)
          }
          
          if (relative_freq_priors) {
            
            clip_db_sub.hapBy <- clip_db_sub.G[clip_db_sub.G[, toHap_col_tmp] %in% rownames(tmp), ]
            toHap_priors <- table(clip_db_sub.hapBy[, toHap_col_tmp])/sum(table(clip_db_sub.hapBy[, toHap_col_tmp]))
            
            if (length(toHap_priors) != nrow(tmp)) {
              toHap_priors_tmp <- c(rep(0, nrow(tmp)))
              names(toHap_priors_tmp) <- rownames(tmp)
              for (i in names(toHap_priors)) {
                toHap_priors_tmp[i] <- toHap_priors[i]
              }
              
              toHap_priors <- toHap_priors_tmp
            }
            
            hap.df <- createHaplotypeTable(tmp, HapByPriors = hapBy_priors, toHapByCol = TRUE, toHapPriors = toHap_priors)
            
            relFreqDf.tmp <- data.frame(c(sample_name, G, hap.df[, 2:length(hap.df)]),
                                        stringsAsFactors = F)
            relFreqDf.tmp <- asNum(relFreqDf.tmp)
            return(relFreqDf.tmp)
          } else {
            hap.df <- createHaplotypeTable(tmp)
            relFreqDf.tmp <- data.frame(c(sample_name, G, hap.df[, 2:length(hap.df)]),
                                        stringsAsFactors = F)
            relFreqDf.tmp <- asNum(relFreqDf.tmp)
            return(relFreqDf.tmp)
            
          }
        }
      }
    }
    
    
    ),use.names=FALSE) %>% as.data.frame()

    colnames(GENES.df.num) <- c("SUBJECT", "GENE", gsub(pattern = "*", "_", hapBy_alleles, fixed = T),
                                'ALLELES', 'PRIORS_ROW', 'PRIORS_COL','COUNTS1', 'K1',
                                'COUNTS2', 'K2','COUNTS3', 'K3','COUNTS4', 'K4')
    
    # Check if toHap_col genes are in toHap_GERM
    if (length(GENES.df.num) == 0)
      stop("Genes in haplotype column to be infered do not match the genes germline given")
    
    
    ## Fill deleted according to a k thershold
    unkIDX <- which(GENES.df.num[gsub("*", "_", hapBy_alleles[1], fixed = T)][, 1] == "Unk")
    delIDX <- unkIDX[which(sapply(unkIDX, function(i) min(GENES.df.num[i,paste0('K',1:4)], na.rm = T)) >= kThreshDel)]
    GENES.df.num[delIDX, paste(gsub("*", "_", hapBy_alleles[1], fixed = T))] <- "Del"
    
    unkIDX <- which(GENES.df.num[gsub("*", "_", hapBy_alleles[2], fixed = T)][, 1] == "Unk")
    delIDX <- unkIDX[which(sapply(unkIDX, function(i) min(GENES.df.num[i,paste0('K',1:4)], na.rm = T)) >= kThreshDel)]
    GENES.df.num[delIDX, paste(gsub("*", "_", hapBy_alleles[2], fixed = T))] <- "Del"
    
    ## Add as unknown genes that do not appear in the individual and mark them as unknown
    
    GENES.MISSING <- GENES.ref[!(GENES.ref %in% GENES.df.num$GENE)]
    if (length(GENES.MISSING) > 0) {
      m <- length(GENES.MISSING)
      
      sub.df <- data.frame(do.call(rbind, lapply(1:m, function(x) {
        c(sample_name, NA, "Unk", "Unk", rep(NA, 11))
      })))
      names(sub.df) <- names(GENES.df.num)
      
      sub.df$GENE <- GENES.MISSING
      
      sub.df[,gsub("*", "_", hapBy_alleles, fixed = T)] <- matrix(rep(ifelse(GENES.MISSING %in% nonReliable_Vgenes_vec,'NR',ifelse(GENES.MISSING %in% deleted_genes_vec, 'Del', 'Unk')), 2), ncol=2)
      
      GENES.df.num <- rbind(GENES.df.num, sub.df)
    }
    
    haplo_db <- rbind(haplo_db, GENES.df.num)
  }
  
  return(haplo_db)
  
}
