# RAbHIT functions -----------------------------------------------------

#' @include rabhit.R
#' @include internal_functions.R
NULL

##########################################################################
#' Anchor gene haplotype inference
#'
#' The \code{createFullHaplotype} functions infers haplotype based on an anchor gene.
#'
#'
#' @param    clip_db               a \code{data.frame} in AIRR format. See details.
#' @param    toHap_col             a vector of column names for which a haplotype should be inferred. Default is v_call and d_call
#' @param    hapBy_col             column name of the anchor gene. Default is j_call
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
#'  \item \code{subject}:        the subject name.
#'  \item \code{gene}:           the gene name.
#'  \item Anchor gene allele 1:  the haplotype inference for chromosome one. The column name is the anchor gene with the first allele.
#'  \item Anchor gene allele 2:  the haplotype inference for chromosome two. The column name is the anchor gene with the second allele.
#'  \item \code{alleles}:        allele calls for the gene.
#'  \item \code{proirs_row}:     priors based on relative allele usage of the anchor gene.
#'  \item \code{proirs_col}:     priors based on relative allele usage of the inferred gene.
#'  \item \code{counts1}:        the appereance count on each chromosome of the first allele from \code{alleles}, the counts are seperated by a comma.
#'  \item \code{k1}:             the Bayesian factor value for the first allele (from \code{alleles}) inference.
#'  \item \code{counts2}:        the appereance count on each chromosome of the second allele from \code{alleles}, the counts are seperated by a comma.
#'  \item \code{k2}:             the Bayesian factor value for the second allele (from \code{alleles}) inference.
#'  \item \code{counts3}:        the appereance count on each chromosome of the third allele from \code{alleles}, the counts are seperated by a comma.
#'  \item \code{k3}:             the Bayesian factor value for the third allele (from \code{alleles}) inference.
#'  \item \code{counts4}:        the appereance count on each chromosome of the fourth allele from \code{alleles}, the counts are seperated by a comma.
#'  \item \code{k4}:             the Bayesian factor value for the fourth allele (from \code{alleles}) inference.
#'}
#'
#' @details
#' Function accepts a \code{data.frame} in AIRR format (\url{https://changeo.readthedocs.io/en/stable/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{'subject'}: The subject name
#'   \item \code{'v_call'}: V allele call(s) (in an IMGT format)
#'   \item \code{'d_call'}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{'j_call'}: J allele call(s) (in an IMGT format)
#' }
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db, HVGERM, HDGERM)
#'
#' # Selecting a single individual
#' clip_db = samples_db[samples_db$subject=='I5', ]
#'
#' # Infering haplotype
#' haplo_db = createFullHaplotype(clip_db,toHap_col=c('v_call','d_call'),
#' hapBy_col='j_call',hapBy='IGHJ6',toHap_GERM=c(HVGERM,HDGERM))
#'
#'
#' @export

createFullHaplotype <- function(clip_db, toHap_col = c("v_call", "d_call"), hapBy_col = "j_call", hapBy = "IGHJ6", toHap_GERM, relative_freq_priors = TRUE,
                                kThreshDel = 3, rmPseudo = TRUE, deleted_genes = c(), nonReliable_Vgenes = c(), min_minor_fraction = 0.3, chain = c("IGH", "IGK", "IGL", "TRB")) {
  
  
  # Check if germline was inputed
  if (missing(toHap_GERM))
    stop("Missing toHap_GERM, please input germline sequences")
  
  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)
  
  if (!("subject" %in% names(clip_db))) {
    
    clip_db$subject <- "S1"
  }
  
  haplo_db <- c()
  clip_db <- clip_db %>% select(.data$subject, !!eval(c(hapBy_col,toHap_col)))
  for (sample_name in unique(clip_db$subject)) {
    
    if (is.list(nonReliable_Vgenes)){
      nonReliable_Vgenes_vec <- nonReliable_Vgenes[[sample_name]]
    } else nonReliable_Vgenes_vec <- nonReliable_Vgenes
    
    if (is.data.frame(deleted_genes)) {
      deleted_genes_vec <- deleted_genes %>% filter(.data$subject == sample_name, .data$deletion == "Deletion") %>% select(.data$gene) %>% pull()
      if (is.null(nonReliable_Vgenes_vec))
        nonReliable_Vgenes_vec <- deleted_genes %>% filter(.data$subject == sample_name, .data$deletion == "Non reliable") %>% select(.data$gene) %>% pull()
    } else deleted_genes_vec <- c()
    
    
    ### Check if haplotype can be infered by the specific gene in the data set.  Only relevant genes with one assignment
    clip_db_sub <- clip_db %>% filter(.data$subject==sample_name, !grepl(',', !!as.name(hapBy_col), perl = T))
    hapBy_priors <- clip_db_sub %>% filter(grepl(paste0(hapBy, "\\*"), !!as.name(hapBy_col), perl = T)) %>% dplyr::count(!!as.name(hapBy_col)) %>% mutate(freq = n / sum(n)) %>% select(-n)
    hapBy_priors <- setNames(hapBy_priors[[2]], hapBy_priors[[1]])
    hapBy_alleles <- names(hapBy_priors)
    if (length(hapBy_alleles) != 2){
      if(sample_name == tail(unique(clip_db$subject),1)){
        stop("Can not haplotype by more or less than two alleles")
      }else{
        message(paste0("For sample ",sample_name,", there were ", length(hapBy_alleles), " alleles, can not haplotype by ", ifelse(length(hapBy_alleles)>2, "more","less")," than two alleles."))
        next()}
    }else{
      message(paste0("For sample ",sample_name, ", haplotyping with ", paste0(hapBy_alleles, collapse = "/")))
    }
    
    if(min(hapBy_priors) < min_minor_fraction){
      if(sample_name == tail(unique(clip_db$subject),1)){
        stop("Can not haplotype, minor allele fraction lower than the cutoff set by the user")
      }else{
        message(paste0("minor allele fraction lower than the cutoff set by the user for sample ",sample_name,", try changing the parameters"))
        next()
      }
    }
    
    
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
        
        toHap_col_tmp <-   toHap_col[stringi::stri_detect_fixed(pattern = substr(tolower(G),4,4),str = toHap_col)]
        
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
    
    colnames(GENES.df.num) <- c("subject", "gene", gsub(pattern = "*", "_", hapBy_alleles, fixed = T),
                                'alleles', 'proirs_row', 'proirs_col','counts1', 'k1',
                                'counts2', 'k2','counts3', 'k3','counts4', 'k4')
    # Check if toHap_col genes are in toHap_GERM
    if (length(GENES.df.num) == 0)
      stop("Genes in haplotype column to be infered do not match the genes germline given")
    
    
    ## Fill deleted according to a k thershold
    unkIDX <- which(GENES.df.num[gsub("*", "_", hapBy_alleles[1], fixed = T)][, 1] == "Unk")
    delIDX <- unkIDX[which(sapply(unkIDX, function(i) min(GENES.df.num[i,paste0('k',1:4)], na.rm = T)) >= kThreshDel)]
    GENES.df.num[delIDX, paste(gsub("*", "_", hapBy_alleles[1], fixed = T))] <- "Del"
    
    unkIDX <- which(GENES.df.num[gsub("*", "_", hapBy_alleles[2], fixed = T)][, 1] == "Unk")
    delIDX <- unkIDX[which(sapply(unkIDX, function(i) min(GENES.df.num[i,paste0('k',1:4)], na.rm = T)) >= kThreshDel)]
    GENES.df.num[delIDX, paste(gsub("*", "_", hapBy_alleles[2], fixed = T))] <- "Del"
    
    ## Add as unknown genes that do not appear in the individual and mark them as unknown
    
    GENES.MISSING <- GENES.ref[!(GENES.ref %in% GENES.df.num$gene)]
    if (length(GENES.MISSING) > 0) {
      m <- length(GENES.MISSING)
      
      sub.df <- data.frame(do.call(rbind, lapply(1:m, function(x) {
        c(sample_name, NA, "Unk", "Unk", rep(NA, 11))
      })))
      names(sub.df) <- names(GENES.df.num)
      
      sub.df$gene <- GENES.MISSING
      
      sub.df[,gsub("*", "_", hapBy_alleles, fixed = T)] <- matrix(rep(ifelse(GENES.MISSING %in% nonReliable_Vgenes_vec,'NR',ifelse(GENES.MISSING %in% deleted_genes_vec, 'Del', 'Unk')), 2), ncol=2)
      
      GENES.df.num <- rbind(GENES.df.num, sub.df)
    }
    
    haplo_db <- rbind(haplo_db, GENES.df.num)
  }
  
  return(haplo_db)
  
}
