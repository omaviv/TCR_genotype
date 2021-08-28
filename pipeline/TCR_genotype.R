library(seqinr)
library(tigger)
library(optparse)
library(gtools)
library(ggplot2)
library(stringr)
library(dplyr)
library(data.table)
library(alakazam)
library(reshape2)
library(utils)
library(rabhit)
library(cowplot)
library(fastmatch)
suppressMessages(library(parallel))
suppressMessages(library(foreach))

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
script.dirname <- ifelse(substr(script.dirname,1,1)=="/", script.dirname, paste0(getwd(), "/",script.dirname))


load(paste0(script.dirname,"/sysdata.rda"))
source(paste0(script.dirname,"/haplotype/createFullHaplotype.R"))
source(paste0(script.dirname,"/haplotype/internal_functions.R"))
# source(paste0(script.dirname,"TIgGER/internal_functions.R"))

####################################################################################################################
######################### PARAMETERS TO INTIAL BEFORE FIRST RUN ####################################################
####################################################################################################################

# The number of threads to use for IgBlast.
num_of_threads <- 44

# max V position to look for novel SNPs 
max_snp_position <- 316

# undocumented_alleles_2_ignore <- c()
undocumented_alleles_2_ignore <- c("TRBV13*01_A170T", "TRBV13*01_T158C", "TRBV10-3*02_C225G", "TRBV20-1*01_C142A", "TRBV30*01_A113C", "TRBV6-6*01_C261T",
                                   "TRBV7-9*05_A19G_C256T",
                                   # "TRBV15*bp02_A316C", "TRBV5-4*bp01_C159T", "TRBV6-6*bp03_G216C", "TRBV6-6*bp03_T201C_A202C_G216C", "TRBV6-6*bp03_T231C_C261T",
                                   # "TRBV15*bp02_G153T", "TRBV19*bp01_T310C_G311C_C314T", "TRBV5-4*bp01_G205A", "TRBV5-5*bp01_G232A", 
                                   "TRBV20-1*ap02_T310G", "TRBV7-8*ap01_T295C", "TRBV7-4*ap01_G291C_A297G", "TRBV7-4*ap01_G291C_A297G_C314T", "TRBV7-9*ap01_G313T")

######################## Environment & tools 
# In case of necessity for loading IgBlast environment before run IgBlast, 
# you should write the command of loading the environment at the beggining with ';' follow after it.
#For example:
load_ig_environment <- "export IGDATA=/home/bcrlab/eitan/data/igblast/; LD_LIBRARY_PATH=/private/anaconda3/lib nice -19"

# IgBlast path. In case that igblast can run simply by igblastn command it could be just "igblastn"
igblastn_path <- "/private/tools/igblast/ncbi-igblast-1.16.0/bin/igblastn"

# Igdiscover environment
igdiscover_env <- "source activate /home/bcrlab/peresay/.conda/envs/igdiscover"

# the path for makeblastdb tool
makeblastdb_path <- "/private/tools/makeblastdb"

# MAKEDB comand (nice -19 for lowering the priority)
makedb_com <- "nice -19 MakeDb.py"
makedb_com <- paste0(makedb_com, " igblast")

################################ Germline refernces #########################################################
# V refernces
full_germline_db_V <- paste0(script.dirname,"/igblast_ref/TRBV_ref")
full_TRBV_imgt_ref <- paste0(script.dirname,"/fasta_references/TRBV.fasta")
full_TRBV_ref <- paste0(script.dirname,"/fasta_references/TRBV_ref.fasta")

sc_germline_db_V <- paste0(script.dirname,"/igblast_ref/TRV_ref")
sc_TRBV_imgt_ref <- paste0(script.dirname,"/fasta_references/TRV.fasta")
sc_TRBV_ref <- paste0(script.dirname,"/fasta_references/TRV_ref.fasta")

biomed_germline_db_V <- paste0(script.dirname,"/igblast_ref/BIOMED2_TRBV_ref")
biomed_TRBV_imgt_ref <- paste0(script.dirname,"/fasta_references/BIOMED2_TRBV.fasta")
biomed_TRBV_ref <- paste0(script.dirname,"/fasta_references/BIOMED2_TRBV_ref.fasta")

adaptive_germline_db_V <- paste0(script.dirname,"/igblast_ref/Adaptive_TRBV_ref")
adaptive_TRBV_imgt_ref <- paste0(script.dirname,"/fasta_references/Adaptive_TRBV.fasta")
adaptive_TRBV_ref <- paste0(script.dirname,"/fasta_references/Adaptive_TRBV_ref.fasta")

# D refernces
germline_db_D <- paste0(script.dirname,"/igblast_ref/TRBD_ref")
TRBD_imgt_ref <- paste0(script.dirname,"/fasta_references/TRBD.fasta")
TRBD_ref <- paste0(script.dirname,"/fasta_references/TRBD_ref.fasta")

# J refernces
full_germline_db_J <- paste0(script.dirname,"/igblast_ref/TRBJ_ref")
full_TRBJ_imgt_ref <- paste0(script.dirname,"/fasta_references/TRBJ.fasta")
full_TRBJ_ref <- paste0(script.dirname,"/fasta_references/TRBJ_ref.fasta")

sc_germline_db_J <- paste0(script.dirname,"/igblast_ref/TRJ_ref")
sc_TRBJ_imgt_ref <- paste0(script.dirname,"/fasta_references/TRJ.fasta")
sc_TRBJ_ref <- paste0(script.dirname,"/fasta_references/TRJ_ref.fasta")

biomed_germline_db_J <- paste0(script.dirname,"/igblast_ref/TRBJ_ref")
biomed_TRBJ_imgt_ref <- paste0(script.dirname,"/fasta_references/TRBJ.fasta")
biomed_TRBJ_ref <- paste0(script.dirname,"/fasta_references/TRBJ_ref.fasta")

adaptive_germline_db_J <- paste0(script.dirname,"/igblast_ref/Adaptive_TRBJ_ref")
adaptive_TRBJ_imgt_ref <- paste0(script.dirname,"/fasta_references/Adaptive_TRBJ.fasta")
adaptive_TRBJ_ref <- paste0(script.dirname,"/fasta_references/Adaptive_TRBJ_ref.fasta")

# Gene usage data frame for deletion detection
full_gene_usages <- read.csv(paste0(script.dirname,"/gene_usage/FULL_TRBV_GENE_MIN_USAGE.tab"), sep = "\t", stringsAsFactors = F)
sc_gene_usages <- full_gene_usages
biomed_gene_usages <- read.csv(paste0(script.dirname,"/gene_usage/BIOMED2_TRBV_GENE_MIN_USAGE.tab"), sep = "\t", stringsAsFactors = F)
adaptive_gene_usages <- read.csv(paste0(script.dirname,"/gene_usage/Adaptive_TRBV_GENE_MIN_USAGE.tab"), sep = "\t", stringsAsFactors = F)



################
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="makedb file name", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="reference file name", metavar="character"),
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="output path folder", metavar="character"),
  make_option(c("--filter_chimera"), type="logical", default=F, action = "store_true",
              help="Ignore novel allele candidates with a potential to result by a chimera sequence.", metavar="character"),
  make_option(c("-t", "--sequence_type"), type="character", default="full", 
              help="The sequence type \n
              - full - full length sequencing. \n
              - sc - single cell sequencing. \n
              - bio - sequenced by Biomed-2 primers \n
              - adapt - sequenced by Adaptive Biotechnologies primers", metavar="character"),
  make_option(c("-c", "--constcount"), type="integer", default=1, 
              help="The minimal 'consensus_count' value to filter by for the genotype inferring.",metavar="number"))


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("input file must be supplied", call.=FALSE)
}

if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("sample must be supplied", call.=FALSE)
}


if (!(tolower(opt$sequence_type) %in% c("full", "sc", "bio", "adapt"))){
  print_help(opt_parser)
  stop("sequence_type must be (full / sc / bio / adapt)", call.=FALSE)
}


igblast_input <- opt$file
SAMP <- opt$sample
output_path <- opt$path
sequencing_length <- tolower(opt$sequence_type)
filter_chimera <- opt$filter_chimera
min_consensus_count <- opt$constcount

####################################################################################################################
##################################### PARAMETERS ###################################################################
####################################################################################################################

# const parameters
sample_path <- paste0(output_path, "/", SAMP, "/")
genotypes_path <- paste0(output_path, "/genotypes/")

igblast_constant_args <- paste("-ig_seqtype", "TCR", "-D_penalty", "-1", "-auxiliary_data", "optional_file/human_gl.aux",
                               "-outfmt", "'7 std qseq sseq btop'", "-domain_system", "imgt", "-num_threads ", num_of_threads)

#### MAKEDB PARAMETERS
makedb_com <- "nice -19 MakeDb.py igblast"
makedb_out_dir <- paste0(output_path, "/MAKEDB/")
makedb_additional_parameters <- paste("--format airr", "--extended")

tcrb_repo <- paste0(output_path, "/TCRB_repo.fasta") # The reference sequences into the "REPO" file for MakeDB.


### PERSONAL NOVEL IGBLAST
novel_igblast_input <- paste0(sample_path, SAMP, ".fasta")
personal_novel_fasta <- paste0(sample_path, "with_novel_V_personal_ref.fasta")
personal_novel_germline_db_V <- paste0(sample_path, "/personal_ref/novel_V")

### PERSONAL GENOTYPE PARAMETERS
personal_ref_path <- paste0(sample_path, "personal_ref/")
v_personal_ref <- paste0(personal_ref_path, "V_personal_ref.fasta")
d_personal_ref <- paste0(personal_ref_path, "D_personal_ref.fasta")
j_personal_ref <- paste0(personal_ref_path, "J_personal_ref.fasta")
personal_repo <- paste0(sample_path, "personal_repo.fasta")
personal_fasta <- paste0(sample_path, SAMP, ".fasta")

personal_db_V <- paste0(personal_ref_path, "genotyped_BV")
personal_db_D <- paste0(personal_ref_path, "genotyped_BD")
personal_db_J <- paste0(personal_ref_path, "genotyped_BJ")

personal_igblast_input <- paste0(sample_path, SAMP, ".fasta")
personal_igblast_output <- paste0(sample_path, SAMP, "_genotyped.fmt7")


if (sequencing_length == "full") {
  germline_db_V <- full_germline_db_V
  germline_db_J <- full_germline_db_J
  TRBV_imgt_ref <- full_TRBV_imgt_ref
  TRBJ_imgt_ref <- full_TRBJ_imgt_ref
  TRBV_ref <- full_TRBV_ref
  TRBJ_ref <- full_TRBJ_ref
  gene_usages <- full_gene_usages
} else if (sequencing_length == "sc") {
  germline_db_V <- sc_germline_db_V
  germline_db_J <- sc_germline_db_J
  TRBV_imgt_ref <- sc_TRBV_imgt_ref
  TRBJ_imgt_ref <- sc_TRBJ_imgt_ref
  TRBV_ref <- sc_TRBV_ref
  TRBJ_ref <- sc_TRBJ_ref
  gene_usages <- sc_gene_usages
} else if (sequencing_length == "bio"){
  germline_db_V <- biomed_germline_db_V
  germline_db_J <- biomed_germline_db_J
  TRBV_imgt_ref <- biomed_TRBV_imgt_ref
  TRBJ_imgt_ref <- biomed_TRBJ_imgt_ref
  TRBV_ref <- biomed_TRBV_ref
  TRBJ_ref <- biomed_TRBJ_ref
  gene_usages <- biomed_gene_usages
} else if (sequencing_length == "adapt"){
  germline_db_V <- adaptive_germline_db_V
  germline_db_J <- adaptive_germline_db_J
  TRBV_imgt_ref <- adaptive_TRBV_imgt_ref
  TRBJ_imgt_ref <- adaptive_TRBJ_imgt_ref
  TRBV_ref <- adaptive_TRBV_ref
  TRBJ_ref <- adaptive_TRBJ_ref
  gene_usages <- adaptive_gene_usages
} else {
  print(paste0("Error 'sequencing_length' value: ", sequencing_length))
  quit()
}
################### Create essential directories for the output files ################
setwd(output_path)
dir.create(SAMP)
dir.create(paste0(output_path, "/igblast/"))
dir.create(makedb_out_dir)
dir.create(sample_path, 'novel_ref')
dir.create(genotypes_path)
# Create the personal reference directory
dir.create(personal_ref_path)

setwd(SAMP)

####################################################################################################################
###################################### READ IMGT GERMLINE REFERENCE ################################################
####################################################################################################################

TRBV_GERM <- read.fasta(TRBV_imgt_ref, as.string = TRUE)
if (grepl("|", names(TRBV_GERM), fixed = T)[1]){
  names(TRBV_GERM) <- sapply(strsplit(names(TRBV_GERM),"|",fixed = T),"[",2) 
}
TRBV_GERM <- toupper(unlist(TRBV_GERM))

TRBD_GERM <- read.fasta(TRBD_imgt_ref, as.string = TRUE)
names(TRBD_GERM) <- sapply(strsplit(names(TRBD_GERM),"|",fixed = T),"[",2) 
TRBD_GERM <- toupper(unlist(TRBD_GERM))

TRBJ_GERM <- read.fasta(TRBJ_imgt_ref, as.string = TRUE)
if (grepl("|", names(TRBJ_GERM), fixed = T)[1]){
  names(TRBJ_GERM) <- sapply(strsplit(names(TRBJ_GERM),"|",fixed = T),"[",2) 
}
TRBJ_GERM <- toupper(unlist(TRBJ_GERM))
# write the reference sequences into the "REPO" file for MakeDB.
write.fasta(sequences = as.list(c(TRBV_GERM, TRBJ_GERM)), names = c(names(TRBV_GERM), names(TRBJ_GERM)),
            tcrb_repo, open="w")


####################################################################################################################
######################################## FIRST IGBLAST + MAKEDB RUNNING ############################################
####################################################################################################################
print("FIRST IGBLAST RUNNING")

pathname = paste0(output_path,"/", SAMP)
igblast_output <- paste0(output_path, "/igblast/",SAMP,".fmt7")
makedb_output <- paste0(makedb_out_dir, SAMP, ".tab")

system(paste(load_ig_environment, igblastn_path,"-germline_db_V", germline_db_V,
             "-germline_db_D", germline_db_D, "-germline_db_J", germline_db_J, igblast_constant_args,
             "-query", igblast_input, "-out", igblast_output,
             "-num_alignments_V 1", "-num_alignments_D 1", "-num_alignments_J 1"))

print("RUNNING MAKEDB ON THE FIRST IGBLAST OUTPUT")
system(paste(makedb_com, "-i", igblast_output, "-s", igblast_input, '-r', tcrb_repo, "-o", makedb_output))

############################################## READ DATA ####################################################################

DATA <- read.delim(makedb_output, header=T, sep = "\t", stringsAsFactors = F)

DATA$v_start <- stringi::stri_locate(DATA$sequence_alignment,regex = "[ATCG]")

allele_diff <- function(germs){
  germs <- lapply(germs, function(x) strsplit(x,'')[[1]])
  germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))
  setdiff_mat <- function(x){
    sum(!unique(x) %in% c('.',NA,"N"))#
  }
  idx = which(apply(germs_m,2,setdiff_mat)>1)
  return(idx)
}

# filter by 3 mutations over the V
v_seqs <- sapply(1:nrow(DATA),function(x) substr(DATA$sequence_alignment[x],1,DATA$v_germline_end[x])) # sequences by length of V gene length found
#Export cluster
cluster <- parallel::makeCluster(15, type="PSOCK")
parallel::clusterExport(cl=cluster,c("allele_diff","DATA","v_seqs","TRBV_GERM"), envir=environment())
DATA$v_mut <- parSapply(cluster,1:nrow(DATA), function(i){
  allele <- DATA$v_call[i];
  allele <- sapply(strsplit(allele, ",", fixed = T), "[", 1);
  idx <- allele_diff(c(v_seqs[i], TRBV_GERM[[allele]]));
  v_min <- min(DATA$v_start[grep(allele, DATA$v_call,fixed=T)])+5;
  sum(idx>v_min & idx<=316)<=3;
}) # get mutations
stopCluster(cluster)

DATA <- DATA[DATA$v_mut,]

####################################################################################################################
####################################### SEARCHING FOR NOVEL ALLELLES ###############################################
####################################################################################################################

print("SEARCHING FOR NOVEL ALLELLES")

convert_format <- function(igd,vgerm){
  #if(grepl("del|ins",igd)){ return(NA)}
  if(igd==""){ return(NA)}
  tmp <- unlist(strsplit(igd,'; '))
  tmp <- tmp[!grepl("del|ins",tmp)]
  mut<- sapply(tmp,function(x){i <- gregexpr("[0-9]", x)[[1]][length(gregexpr("[0-9]", x)[[1]])];
  m <- as.numeric(substr(x,start = 1,stop = i) );
  vgerm.nogap <- vgerm[vgerm!='.']  ;
  vgerm.nogap.str <- substr(paste0(vgerm.nogap,collapse = ""),1,m);
  toadd <- sum(s2c(sub("[.]*$", "", togap(paste0(vgerm,collapse = ""),vgerm.nogap.str),fixed = F,perl = T))=='.')
  m <- m + toadd;
  n <- unlist(strsplit('>',x=substr(x,start = i+1,stop = nchar(x)),fixed=T));
  paste0(n[1],m,n[2])})
  paste(mut,collapse = '_')
}

togap <- function(vgap,vdj){                      ##add in vdj gaps
  gapadd <- vdj
  for(i in which(unlist(strsplit(vgap,"",fixed=T)) == ".") ){
    gapadd <- paste0(substr(gapadd,1,i-1),".",substr(gapadd,i,nchar(gapadd)))
  }
  return(gapadd)
}


pathname = paste0(output_path,"/", SAMP)
dir.create(paste0(output_path,"/", SAMP,"/database"))

file.copy(TRBV_ref, paste0(output_path,"/", SAMP,"/database/V.fasta"), overwrite = T)
file.copy(TRBD_ref, paste0(output_path,"/", SAMP,"/database/D.fasta"), overwrite = T)
file.copy(TRBJ_ref, paste0(output_path,"/", SAMP,"/database/J.fasta"), overwrite = T)


system(paste0(igdiscover_env," && cd ",pathname,"; ",
              'igdiscover igblast --threads 1 --no-cache',
              ' --pathig ', igblast_output,
              ' database ', igblast_input,
              " > ", SAMP, "_igdiscover-pass.tab"))

data_igdiscover <- read.delim(paste0(SAMP, "_igdiscover-pass.tab"), sep = '\t', header = T, stringsAsFactors = F)
data_igdiscover$sequence_id <- sapply(data_igdiscover$name, function(x) strsplit(x,"|",fixed=T)[[1]][1])
data_igdiscover <- data_igdiscover[data_igdiscover$sequence_id %in% DATA$sequence_id,]

ig_ind <- 1
for (seq_id in data_igdiscover$sequence_id) {
  seq_ind <- which(DATA$sequence_id==seq_id)
  cdr3 <- DATA$junction[[seq_ind]]
  if (str_length(cdr3) > 6) {
    cdr3 <- substr(cdr3, 4, str_length(cdr3)-3)
    cdr3_aa <- DATA$junction_aa[[seq_ind]]
    cdr3_aa <- substr(cdr3_aa, 2, str_length(cdr3_aa)-1)
  } else {
    cdr3 <- ""
    cdr3_aa <- ""
  }
  data_igdiscover$CDR3_nt[ig_ind] <- cdr3
  data_igdiscover$CDR3_aa[ig_ind] <- cdr3_aa
  ig_ind <- ig_ind + 1
}


write.table(data_igdiscover, paste0(SAMP, "_igdiscover-pass_mut.tab"), sep = "\t")

system(paste0(igdiscover_env," && cd ",pathname,"; ",
              paste0('igdiscover discover --threads ', num_of_threads,
                     ' --consensus-threshold 50 --database database/V.fasta ',
                     ' -o ./ ',SAMP,'_igdiscover-pass_mut.tab',
                     " > ", SAMP, "_candidates_igdiscover.tab")))

novel_igdiscover <- read.delim(paste0(SAMP, "_candidates_igdiscover.tab"), sep = '\t', header = T, stringsAsFactors = F)

novel_igdiscover <- novel_igdiscover[novel_igdiscover$database_diff!=0,]
TRBV_GERM_NEW <- TRBV_GERM
if(nrow(novel_igdiscover)>0) novel_igdiscover$NT_SUBSTITUTIONS <- sapply(1:nrow(novel_igdiscover),
                                                                         function(i){
                                                                           convert_format(novel_igdiscover$database_changes[i],
                                                                                          vgerm = seqinr::s2c(TRBV_GERM_NEW[[novel_igdiscover$source[i]]]))})

if(nrow(novel_igdiscover)>0)  novel_igdiscover <- novel_igdiscover[!is.na(novel_igdiscover$NT_SUBSTITUTIONS),]
if(nrow(novel_igdiscover)>0){
  
  novel_igdiscover <- novel_igdiscover[grep('_',novel_igdiscover$name,fixed = T),]
  
  ### filter inferences that have at least 2 CDR3 length and at least 2 TRBJ genes
  novel_igdiscover <- novel_igdiscover[novel_igdiscover$CDR3s>=2 & novel_igdiscover$Js>=2,]
}
if(nrow(novel_igdiscover)>0){
  ### filter positions out of range
  # min range is N+5 where N is the position of the primer end
  # max range is 316 for all
  novel_igdiscover$MIN_V_START <- sapply(1:nrow(novel_igdiscover),function(i) min(DATA$v_start[grep(novel_igdiscover$source[i], DATA$v_call,fixed=T)])+5)
  novel_igdiscover$NT_SUBSTITUTIONS_OR <-  novel_igdiscover$NT_SUBSTITUTIONS
  novel_igdiscover$NT_SUBSTITUTIONS <- sapply(1:nrow(novel_igdiscover), function(i){
    
    snps <- novel_igdiscover$NT_SUBSTITUTIONS[i]
    subs <- strsplit(snps,'_')[[1]]
    allele <- novel_igdiscover$source[i]
    
    # filter minimum
    min_start_pos <- novel_igdiscover$MIN_V_START[i]
    positions <- sapply(subs,str_extract,pattern = "[0-9]+")
    positions <- as.integer(positions)
    idx_min <- which(positions < min_start_pos)
    
    
    # filter max
    idx_max <- which(positions > max_snp_position)
    
    
    idx <- c(idx_min,idx_max)
    if(length(idx)!=0) subs <- subs[-idx]
    
    if(length(subs)!=0) return(paste0(subs,collapse = "_"))
    else return(NA)
  })
  
  novel_igdiscover <- novel_igdiscover[!is.na(novel_igdiscover$NT_SUBSTITUTIONS),]
}
if(nrow(novel_igdiscover)>0){
  novel_igdiscover$POLYMORPHISM_CALL <- paste0(novel_igdiscover$source,"_",novel_igdiscover$NT_SUBSTITUTIONS)
  novel_igdiscover <- novel_igdiscover %>% mutate(GENE = alakazam::getGene(POLYMORPHISM_CALL, strip_d = F))
  
  novel_igdiscover$NOVEL_IMGT <- unlist(sapply(1:nrow(novel_igdiscover),function(i){
    seq <- TRBV_GERM_NEW[[novel_igdiscover$source[i]]]
    snps <- novel_igdiscover$NT_SUBSTITUTIONS[i]
    subs <- strsplit(snps,'_')[[1]]
    positions <- sapply(subs,str_extract,pattern = "[0-9]+")
    positions <- as.integer(positions)
    
    for(ii in 1:length(subs)){
      substr(seq,positions[ii],positions[ii]) <- strsplit(subs[ii],positions[ii])[[1]][2]
    }
    return(seq)
  }))
  
  novel_igdiscover$NOTE <- ""
  
  
  # Remove potential SNPs that close to the 3' end of the V segment or identical sequences between N+5 to 316
  rows2remove <- c()
  max_snp_position <- 316
  for (i in 1:nrow(novel_igdiscover)){
    
    ### if no exact matches for the novel allele then remove.
    ## update exact number to the snps limit
    seq <- substr(novel_igdiscover$NOVEL_IMGT[i],novel_igdiscover$MIN_V_START[i],max_snp_position)
    allele <- strsplit(novel_igdiscover$POLYMORPHISM_CALL[[i]],"_")[[1]][1]
    
    novel_igdiscover$full_exact_new[i] <- sum(grepl(seq, DATA$sequence_alignment, fixed = T)&grepl(allele, DATA$v_call, fixed = T))
    
    if(novel_igdiscover$full_exact_new[i]<(length(grep(allele,DATA$v_call,fixed=T))*0.05)){
      rows2remove <- c(rows2remove, i)
      novel_igdiscover$NOTE[i] <- paste0("Less than 5%(",(length(grep(allele,DATA$v_call,fixed=T))*0.05),") of the gene sequences were exact copies;")
      #next()
    }else{
      novel_igdiscover$NOTE[i] <- paste0("More than 5%(",(length(grep(allele,DATA$v_call,fixed=T))*0.05),") of the gene sequences were exact copies;")
    }
    
    gene <- unlist(str_split(novel_igdiscover$POLYMORPHISM_CALL[[i]], "[*]"))[1]
    gene <- unlist(str_split(gene, "-"))[1] # greps the family?
    ALLELES <- TRBV_GERM[grepl(gene, names(TRBV_GERM))]
    novel_imgt_seq <- novel_igdiscover$NOVEL_IMGT[[i]]
    for (j in 1:length(ALLELES)) {
      if (grepl(ALLELES[[j]], novel_imgt_seq)) {
        new_name <- paste0(names(ALLELES[j]), "_", 
                           gsub(TRBV_GERM[names(ALLELES[j])], 
                                as.character(str_length(TRBV_GERM[names(ALLELES[j])])+1), 
                                novel_imgt_seq), 
                           str_length(novel_imgt_seq))
        name_index <- names(TRBV_GERM) == names(ALLELES[j])
        if (!grepl("_", names(ALLELES[j]))) {
          TRBV_GERM[names(ALLELES[j])] <- novel_imgt_seq
          names(TRBV_GERM)[name_index] <- new_name
        }
        rows2remove <- c(unlist(rows2remove), i)
        novel_igdiscover$NOTE[i] <- paste0(novel_igdiscover$NOTE[i],"Extantion of known allele;")
        next()
      }
      else {
        cutted_allele_seq <- substr(ALLELES[[j]], 1, max_snp_position)
        cutted_allele_seq <- gsub(".", "", cutted_allele_seq, fixed = T)
        cutted_allele_seq <- substr(cutted_allele_seq, 5, str_length(cutted_allele_seq))
        
        cutted_novel <- substr(novel_imgt_seq, 1, max_snp_position)
        cutted_novel <- gsub(".", "", cutted_novel, fixed = T)
        cutted_novel <- substr(cutted_novel, 5, str_length(cutted_novel))
        if(cutted_allele_seq == cutted_novel) {
          rows2remove <- c(unlist(rows2remove), i)
          novel_igdiscover$NOTE[i] <- paste0(novel_igdiscover$NOTE[i],"Identical to known allele;")
        }
      }
      # if (str_length(novel_imgt_seq) > max_snp_position) { 
      #   if (grepl(substr(novel_imgt_seq, 1, max_snp_position), ALLELES[[j]])) {
      #     rows2remove <- c(unlist(rows2remove), i)
      #     novel_igdiscover$NOTE[i] <- paste0(novel_igdiscover$NOTE[i],"Extantion of known allele;")
      #   }
      # }
    }
  }
  
  
  
  if (length(rows2remove)) {
    write.table(novel_igdiscover,file = file.path(pathname,paste0(SAMP,"_pre_novel_selected_igdiscover.tab")), sep='\t', row.names = F)
    novel_igdiscover <- novel_igdiscover[-rows2remove,]
    write.table(novel_igdiscover,file = file.path(pathname,paste0(SAMP,"_novel_selected_igdiscover.tab")), sep='\t', row.names = F)
  }else{
    write.table(novel_igdiscover,file = file.path(pathname,paste0(SAMP,"_novel_selected_igdiscover.tab")), sep='\t', row.names = F)
  }
}

new_novel_df_H <- novel_igdiscover


############################################### filter chimera sequences ###############################################
if (filter_chimera) {
  known_alleles <- TRBV_GERM
  novel_alleles <- new_novel_df_H$NOVEL_IMGT
  novel_alleles <- sapply(novel_alleles,function(x){gsub('-','.',x,fixed = T)})
  names(novel_alleles) <- new_novel_df_H$POLYMORPHISM_CALL
  if (length(novel_alleles) > 0) {
    novel_allele_mismatches <- list()
    for (novel_allele in names(novel_alleles)) {
      if (grepl("_[0-9]", novel_allele)) {
        next
      }
      gene <- sapply(strsplit(novel_allele,"*", fixed = T), '[', 1)
      gene_family <- sapply(strsplit(gene,"-", fixed = T), '[', 1)
      family_alleles <- known_alleles[grepl(gene_family, names(known_alleles))]
      
      novel_allele_mismatches[[novel_allele]] <- list()
      novel_allele_seq <- unlist(strsplit(as.character(novel_alleles[novel_allele]), ""))
      for (known_allele in names(family_alleles)) {
        known_allele_seq <- unlist(strsplit(as.character(family_alleles[known_allele]), ""))
        mismatches_counter <- c(0)
        for (pos in 2:max_snp_position) {
          mismatches_counter[pos] <- mismatches_counter[pos-1]
          if (pos > min(str_length(novel_allele_seq), str_length(known_allele_seq))) {
            next
          }
          if ((novel_allele_seq[[pos]] != ".") & (known_allele_seq[[pos]] != ".") & (novel_allele_seq[[pos]] != known_allele_seq[[pos]])) {
            mismatches_counter[pos] <- mismatches_counter[pos] + 1
          }
        }
        novel_allele_mismatches[[novel_allele]][[known_allele]] <- list()
        novel_allele_mismatches[[novel_allele]][[known_allele]][["prefix"]] <- mismatches_counter
        novel_allele_mismatches[[novel_allele]][[known_allele]][["suffix"]] <- mismatches_counter[max_snp_position] - mismatches_counter
      }
    }
    
    
    chimera_alleles <- c()
    prefix_alleles <- c()
    suffix_alleles <- c()
    min_mismatches <- c()
    position <- c()
    for (novel_allele in names(novel_allele_mismatches)) {
      for (prefix_allele in names(novel_allele_mismatches[[novel_allele]])) {
        for (suffix_allele in names(novel_allele_mismatches[[novel_allele]])) {
          if (suffix_allele != prefix_allele) {
            chimera_alleles <- c(chimera_alleles, novel_allele)
            prefix_alleles <- c(prefix_alleles, prefix_allele)
            suffix_alleles <- c(suffix_alleles, suffix_allele)
            mismatches <- novel_allele_mismatches[[novel_allele]][[prefix_allele]][["prefix"]] + novel_allele_mismatches[[novel_allele]][[suffix_allele]][["suffix"]]
            min_mismatch <- min(mismatches)
            min_mismatches <- c(min_mismatches, min_mismatch)
            position <- c(position, which.min(mismatches))
          }
        }
      }
    }
    
    chimera_df <- data.frame(chimera_alleles, prefix_alleles, suffix_alleles, min_mismatches, position, stringsAsFactors = F)
    if (nrow(chimera_df)) {
      chimera_df$snp_count <- str_count(chimera_df$chimera_alleles, "_")
      chimera_df <- chimera_df[chimera_df$min_mismatches < chimera_df$snp_count,]
    }
    
    if (nrow(chimera_df)) {
      pos_chimera_df <- do.call(rbind, unname(by(chimera_df, chimera_df$chimera_alleles, function(x) x[x$min_mismatches == min(x$min_mismatches),])))
      pos_chimera_df$prefix_gene <- sapply(strsplit(as.character(pos_chimera_df$prefix_alleles), "*", fixed = T), '[', 1)
      pos_chimera_df$suffix_gene <- sapply(strsplit(as.character(pos_chimera_df$suffix_alleles), "*", fixed = T), '[', 1)
      pos_chimera_df <- pos_chimera_df[pos_chimera_df$prefix_gene != pos_chimera_df$suffix_gene,]
      
      prob_chimera <- unique(pos_chimera_df$chimera_alleles[pos_chimera_df$min_mismatches==0])
    }
    else {
      prob_chimera <- list()
    }
  } else {
    prob_chimera <- list()
  }
  
  
  chimeras_df <- data.frame(chimera_full_name=prob_chimera, stringsAsFactors = F)
  if (nrow(chimeras_df)>0) {
    chimeras_df$chimera_seq <- unlist(lapply(chimeras_df$chimera_full_name, function(ch_name){novel_alleles[[ch_name]]}))
    chimeras_df$gene <- sapply(strsplit(chimeras_df$chimera_full_name, "*", fixed = T), "[", 1)
    chimeras_df$new_name <- paste0("ch", 1:nrow(chimeras_df))
    chimeras_df$new_name <- paste(chimeras_df$gene, chimeras_df$new_name, sep = "*")
    novel_alleles <- novel_alleles[!names(novel_alleles) %in% prob_chimera]
    write.table(chimeras_df, file = paste0(sample_path, SAMP, "_chimeras.tsv"),quote = F,row.names = F,sep = "\t")
    
  }
  
  novel_genotype <- novel_alleles
  for (novel_allele in names(novel_genotype)) {
    if (!(novel_allele %in% names(TRBV_GERM))) {
      TRBV_GERM <- c(TRBV_GERM,novel_genotype[novel_allele])
    }
  }
} else {
  ## ADD novel alleles to the genotype database
  novel_genotype <- new_novel_df_H$NOVEL_IMGT
  novel_genotype <- sapply(novel_genotype,function(x){gsub('-','.',x,fixed = T)})
  names(novel_genotype) <- new_novel_df_H$POLYMORPHISM_CALL
  for (novel_allele in names(novel_genotype)) {
    if (!(novel_allele %in% names(TRBV_GERM))) {
      TRBV_GERM <- c(TRBV_GERM,novel_genotype[novel_allele])
    }
  }
}
## Combine the genotyped and others and write to a fasta file for reference
write.fasta(sequences=as.list(gsub(TRBV_GERM, pattern = '.', replacement = '', fixed = T)),
            names=names(TRBV_GERM), personal_novel_fasta, open="w")

## Write the VDJ sequences to realign
print("Collape identical VDJ sequences")

if (!"consensus_count" %in% names(DATA)) {
  if ("reads" %in% names(DATA)) {
    if ("templates" %in% names(DATA)) {
      DATA$templates[!unlist(lapply(DATA$templates, is.integer))] <- 1
      DATA$reads[!unlist(lapply(DATA$reads, is.integer))] <- DATA$templates[!unlist(lapply(DATA$reads, is.integer))]
    } else {
      DATA$reads[!unlist(lapply(DATA$reads, is.integer))] <- 1
      DATA$templates <- 1
    }
    DATA$consensus_count <- DATA$reads
    DATA$duplicate_count <- DATA$templates
  }
  else {
    DATA$consensus_count <- 1
    DATA$duplicate_count <- 1
  }
} else{
  DATA$duplicate_count[!unlist(lapply(DATA$duplicate_count, is.integer))] <- 1
  DATA$consensus_count[!unlist(lapply(DATA$consensus_count, is.integer))] <- DATA$duplicate_count[!unlist(lapply(DATA$consensus_count, is.integer))]
}


DATA <- DATA %>% select(sequence, sequence_alignment, sequence_id,  consensus_count, duplicate_count)
DATA$sequence_alignment <- gsub(".", "", DATA$sequence_alignment, fixed = T)

vdj_seqs <- unique(DATA$sequence_alignment)
for (vdj_seq in vdj_seqs) {
  vdj_indexes <- which(DATA$sequence_alignment==vdj_seq)
  if (length(vdj_indexes) > 1) {
    temp <- DATA[vdj_indexes,]
    temp$consensus_count <- as.numeric(temp$consensus_count)
    temp$duplicate_count <- as.numeric(temp$duplicate_count)
    
    DATA <- DATA[-vdj_indexes,]
    
    seq_id_ind <- which.max(temp$consensus_count);
    seq_id <- temp$sequence_id[seq_id_ind]
    seq_id_ind <- vdj_indexes[seq_id_ind]
    sequence <- temp$sequence[seq_id_ind]
    
    conscount <- sum(temp$consensus_count)
    dupcount <- sum(temp$duplicate_count)
    DATA[nrow(DATA)+1,] <- c(sequence, vdj_seq, seq_id, conscount, dupcount)
  }
}

seq.names <- sapply(1:nrow(DATA),function(x){ paste0(names(DATA)[3:ncol(DATA)],
                                                     rep('=',length(3:ncol(DATA))),
                                                     DATA[x,3:ncol(DATA)],
                                                     collapse = '|')})
seq.names <- gsub('sequence_id=','',seq.names,fixed = T)

write.fasta(sequences=as.list(DATA$sequence), names=seq.names, novel_igblast_input, open="w")
################################## Create a personal reference database file using blast - makeblastdb ###########
pathname <- paste0(output_path, "/", SAMP)
print("CREATING A REFERENCE DATABASE WITH FOUND NOVEL ALLELES")
system(paste(makeblastdb_path,'-parse_seqids -dbtype nucl',
             "-in",  personal_novel_fasta, "-out", personal_novel_germline_db_V))

####################################################################################################################
########################### SECOND IGBLAST + MAKEDB RUNNING WITH PERSONAL NOVEL ALLELES#############################
####################################################################################################################
### Realign using personal reference and IgBlast
print("RUNNING IGBLAST WITH PERSONAL REFERENCE")
pathname = paste0(output_path,"/", SAMP)
makedb_output <- paste0(sample_path, SAMP, "_novel.tab")
novel_repo <- paste0(output_path, "/", SAMP, "/", "novel_repo.fasta")
igblast_output <- paste0(pathname,"/igblast_novel.out")

system(paste(load_ig_environment, igblastn_path, "-germline_db_V", personal_novel_germline_db_V,
             "-germline_db_D", germline_db_D, "-germline_db_J", germline_db_J,
             igblast_constant_args, "-query", novel_igblast_input, "-out", igblast_output))

print("RUNNING MAKEDB ON THE 2ND IGBLAST OUTPUT")
write.fasta(sequences = as.list(c(TRBV_GERM, TRBJ_GERM)), names = c(names(TRBV_GERM), names(TRBJ_GERM)),
            novel_repo, open="w")
system(paste( makedb_com, "-i",  igblast_output, "-s", novel_igblast_input, "-r", novel_repo, "-o", makedb_output))

####################################################################################################################
########################################## GENOTYPE INFERRING ######################################################
####################################################################################################################
DATA <- read.delim(makedb_output, header=T, sep = "\t", stringsAsFactors = F)
# filter by the selected minimal constcount.
DATA <- DATA[DATA$consensus_count >= min_consensus_count,]

# filter by zero mutations over the V
DATA$v_seq <- substr(DATA$sequence_alignment,1,sapply(DATA$v_germline_end, min, max_snp_position)) # sequences by length of V gene length found
DATA$v_mut <- sapply(tigger::getMutCount(DATA$v_seq, DATA$v_call, germline_db = TRBV_GERM),function(x){x[[1]]}) # get mutations
DATA <- DATA[DATA$v_mut <= 1,]

DATA_V_SA <- DATA[!grepl(pattern = ',',DATA$v_call),]
DATA_V_SA <- DATA_V_SA[!DATA_V_SA$v_call %in% undocumented_alleles_2_ignore,]

geno_BV <- inferGenotypeBayesian(DATA_V_SA, germline_db = TRBV_GERM, find_unmutated = FALSE,
                                 novel = new_novel_df_H,v_call = 'v_call')
names(geno_BV) <- names(geno_BV)
geno_BV$GENOTYPED_ALLELES <-apply(geno_BV[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:5]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})

if (filter_chimera) {
  trizygous <- geno_BV[str_count(geno_BV$GENOTYPED_ALLELES, ",") >=2 ,]
  trizygous <- trizygous %>% tidyr::separate_rows(GENOTYPED_ALLELES, sep = ",")
  trizygous$full_name <- paste0(trizygous$gene, "*", trizygous$GENOTYPED_ALLELES)
  novel_names <- trizygous$full_name[grepl("_[A-Z]", trizygous$GENOTYPED_ALLELES)]
  
  if (length(novel_names)){
    novel_allele_mismatches <- list()
    for (novel in novel_names) {
      novel_seq <- unlist(strsplit(as.character(TRBV_GERM[[novel]]), ""))
      gene <- sapply(strsplit(novel, "*", fixed=T), "[", 1)
      gene_alleles <- trizygous$full_name[trizygous$gene == gene]
      gene_alleles <- TRBV_GERM[gene_alleles]
      
      for (allele in names(gene_alleles)) {
        if (allele == novel) {next}
        allele_seq <- unlist(strsplit(as.character(gene_alleles[allele]), ""))
        mismatches_counter <- c(0)
        for (pos in 2:max_snp_position) {
          mismatches_counter[pos] <- mismatches_counter[pos-1]
          if (pos > min(str_length(novel_seq), str_length(allele_seq))) {
            next
          }
          if ((novel_seq[[pos]] != ".") & (allele_seq[[pos]] != ".") & (novel_seq[[pos]] != allele_seq[[pos]])) {
            mismatches_counter[pos] <- mismatches_counter[pos] + 1
          }
        }
        novel_allele_mismatches[[novel]][[allele]] <- list()
        novel_allele_mismatches[[novel]][[allele]][["prefix"]] <- mismatches_counter
        novel_allele_mismatches[[novel]][[allele]][["suffix"]] <- mismatches_counter[max_snp_position] - mismatches_counter
      }
    }
    
    chimera_alleles <- c()
    prefix_alleles <- c()
    suffix_alleles <- c()
    min_mismatches <- c()
    position <- c()
    for (novel_allele in names(novel_allele_mismatches)) {
      for (prefix_allele in names(novel_allele_mismatches[[novel_allele]])) {
        for (suffix_allele in names(novel_allele_mismatches[[novel_allele]])) {
          if (suffix_allele != prefix_allele) {
            chimera_alleles <- c(chimera_alleles, novel_allele)
            prefix_alleles <- c(prefix_alleles, prefix_allele)
            suffix_alleles <- c(suffix_alleles, suffix_allele)
            mismatches <- novel_allele_mismatches[[novel_allele]][[prefix_allele]][["prefix"]] + novel_allele_mismatches[[novel_allele]][[suffix_allele]][["suffix"]]
            min_mismatch <- min(mismatches)
            min_mismatches <- c(min_mismatches, min_mismatch)
            position <- c(position, which.min(mismatches))
          }
        }
      }
    }
    
    chimera_df <- data.frame(chimera_alleles, prefix_alleles, suffix_alleles, min_mismatches, position, stringsAsFactors = F)
    if (nrow(chimera_df)) {
      chimera_df$snp_count <- str_count(chimera_df$chimera_alleles, "_")
      chimera_df <- chimera_df[chimera_df$min_mismatches < chimera_df$snp_count,]
    }
    
    if (nrow(chimera_df)) {
      pos_chimera_df <- do.call(rbind, unname(by(chimera_df, chimera_df$chimera_alleles, function(x) x[x$min_mismatches == min(x$min_mismatches),])))
      pos_chimera_df$prefix_gene <- sapply(strsplit(as.character(pos_chimera_df$prefix_alleles), "*", fixed = T), '[', 1)
      pos_chimera_df$suffix_gene <- sapply(strsplit(as.character(pos_chimera_df$suffix_alleles), "*", fixed = T), '[', 1)
      
      prob_chimera <- unique(pos_chimera_df$chimera_alleles[pos_chimera_df$min_mismatches==0])
    } 
    
    if (length(prob_chimera)) {
      DATA_V_SA <- DATA_V_SA[!DATA_V_SA$v_call %in% prob_chimera,]
      geno_BV <- inferGenotypeBayesian(DATA_V_SA, germline_db = TRBV_GERM, find_unmutated = FALSE,
                                       novel = new_novel_df_H,v_call = 'v_call')
      geno_BV$GENOTYPED_ALLELES <-apply(geno_BV[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:5]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})
    }
  }
}


DATA_D_geno <- DATA[(!grepl(pattern = ',',DATA$d_call) & DATA$d_call !='None') & (DATA$d_sequence_end - DATA$d_sequence_start >= 8),]
DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno$sequence_id),]

# extract d sequence in the direct orientation
DATA_D_reg <- DATA_D_geno[DATA_D_geno$d_germline_start < DATA_D_geno$d_germline_end,]
DATA_D_reg$d_seq <- substr(DATA_D_reg$sequence, DATA_D_reg$d_sequence_start, DATA_D_reg$d_sequence_end)

# extract convert d sequence in the inverted orientation to the direct orientation
DATA_D_inv <- DATA_D_geno[DATA_D_geno$d_germline_start > DATA_D_geno$d_germline_end,]
DATA_D_inv$d_seq <- substr(DATA_D_inv$sequence, DATA_D_inv$d_sequence_start, DATA_D_inv$d_sequence_end)
DATA_D_inv$d_seq <- stringi::stri_reverse(DATA_D_inv$d_seq)
DATA_D_inv$d_seq <- gsub("A", "t", DATA_D_inv$d_seq)
DATA_D_inv$d_seq <- gsub("T", "a", DATA_D_inv$d_seq)
DATA_D_inv$d_seq <- gsub("G", "c", DATA_D_inv$d_seq)
DATA_D_inv$d_seq <- gsub("C", "g", DATA_D_inv$d_seq)
DATA_D_inv$d_seq <- toupper(DATA_D_inv$d_seq)

d_germ_end <- DATA_D_inv$d_germline_start
DATA_D_inv$d_germline_start <- DATA_D_inv$d_germline_end
DATA_D_inv$d_germline_end <- d_germ_end

DATA_D_geno <- rbind(DATA_D_reg, DATA_D_inv)

# filter by zero mutations over the D segment
DATA_D_geno$mut_d <- unlist(lapply(1:nrow(DATA_D_geno), function(i) {
  mut <- 0;
  row_seq <- unlist(strsplit(DATA_D_geno$d_seq[[i]], ""));
  allele_seq <- unlist(strsplit(TRBD_GERM[[DATA_D_geno$d_call[[i]]]], ""));
  for (pos in DATA_D_geno$d_germline_start[[i]]:DATA_D_geno$d_germline_end[[i]]) {
    if (row_seq[pos - (DATA_D_geno$d_germline_start[[i]] - 1)] != allele_seq[pos]) {
      mut <- mut + 1;
    }
  }
  mut
}))

DATA_D_geno <- DATA_D_geno[DATA_D_geno$mut_d == 0,]


geno_BD <- inferGenotypeBayesian(DATA_D_geno, find_unmutated = FALSE, germline_db = TRBD_GERM, v_call = 'd_call')
geno_BD$GENOTYPED_ALLELES <-apply(geno_BD[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})

# handling the D2 alignment error rate, the 01 frequency over heterozygous is between 0.2066 to 0.8969
D2_total <- nrow(DATA_D_geno[grepl("TRBD2", DATA_D_geno$d_call),])
D2_01_count <- nrow(DATA_D_geno[DATA_D_geno$d_call=="TRBD2*01",])
D2_01_freq <- D2_01_count / D2_total

if (D2_01_freq < 0.2066) {
  geno_BD$GENOTYPED_ALLELES[geno_BD$gene == "TRBD2"] <- "02"
} else if (D2_01_freq > 0.8969) {
  geno_BD$GENOTYPED_ALLELES[geno_BD$gene == "TRBD2"] <- "01"
} else if (geno_BD$GENOTYPED_ALLELES[geno_BD$gene == "TRBD2"] == "01") {
  geno_BD$GENOTYPED_ALLELES[geno_BD$gene == "TRBD2"] <- "01,02"
}


DATA_J_SA <- DATA[!grepl(pattern = ',',DATA$j_call),]
geno_BJ <- inferGenotypeBayesian(DATA, germline_db = TRBJ_GERM, find_unmutated = FALSE, v_call = 'j_call')
geno_BJ$GENOTYPED_ALLELES <-apply(geno_BJ[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})

####################################################################################################################
###################################### TRBV DELETION DETECTION #####################################################
####################################################################################################################

DATA <- read.delim(makedb_output, header=T, sep = "\t", stringsAsFactors = F)
# filter by the selected minimal constcount.
DATA <- DATA[DATA$consensus_count >= min_consensus_count,]

cutoff <- 0.0005
p_val_cutoff = 0.001

DATA$v_gene <- unlist(lapply(DATA$v_call, function(Vcall){
  assignments <- unlist(strsplit(Vcall, ",", fixed = T));
  v_genes <- unique(sapply(strsplit(assignments, "*", fixed = T), "[", 1));
  v_genes <- v_genes[order(v_genes)];
  paste(v_genes, collapse = ",")
}))

DATA <- DATA[!grepl(",", DATA$v_gene),]
DATA <- DATA[grepl("TRBV", DATA$v_gene),] # CHeck deletions only for TRBV, we don't have the usage of TRAV yet..

gene_usages$N <- unlist(lapply(gene_usages$GENE, function(gene){nrow(DATA[DATA$v_gene==gene,])}))
gene_usages$TOTAL <- nrow(DATA)
gene_usages$USAGE <- gene_usages$N / gene_usages$TOTAL

gene_usages <- gene_usages[gene_usages$MIN_FREQ != Inf & gene_usages$AVG_USAGE > 1.5*cutoff,]

# calculate the p value for each gene in case that the gene is deleted by binom test
gene_usages$PVAL <- sapply(1:nrow(gene_usages), function(i) {
  if (gene_usages$USAGE[i] < cutoff) {
    return(binom.test(x = gene_usages$N[i], n = gene_usages$TOTAL[i], p = gene_usages$MIN_FREQ[i])$p.value)
  } else {
    return(1)
  }
})

# Detect according to the p values if there are deleted genes
gene_usages$DELETED <- sapply(1:nrow(gene_usages), function(i) {
  if (gene_usages$PVAL[i] <= p_val_cutoff) {
    if ((gene_usages$USAGE[i] < cutoff) & gene_usages$MIN_FREQ[i] != Inf) {
      return(TRUE)
    }
  }
  return(FALSE)
})

gene_usages <- gene_usages[gene_usages$DELETED,]

if (nrow(gene_usages) > 0) {
  deleted_genes <- gene_usages$GENE
  geno_BV <- geno_BV[!geno_BV$gene %in% deleted_genes,]
  
  for (gene in deleted_genes) {
    geno_BV[nrow(geno_BV) + 1,] <- c(gene, NA, NA, NA, NA, NA, NA, NA, NA, 1000, "Deletion")
  }
}


####################################################################################################################
#################################### CREATE PERSONAL GENOTYPE REFERENCE ############################################
####################################################################################################################

## Remove from TRBV_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBV_GERM),'*',fixed=T),'[',1) %in%  geno_BV$gene)
TRBV_GERM.NEW <- TRBV_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_BV)){
  # Remove deleted genes
  if (geno_BV$GENOTYPED_ALLELES[i] == "Deletion") {
    next
  }
  
  gene <- geno_BV$gene[i]
  
  alleles <- geno_BV$GENOTYPED_ALLELES[i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBV_GERM) %in%  paste(gene,alleles,sep='*')
  TRBV_GERM.NEW <- c(TRBV_GERM.NEW,TRBV_GERM[IND])
}


## Remove from TRBD_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBD_GERM),'*',fixed=T),'[',1) %in%  geno_BD$gene)
TRBD_GERM.NEW <- TRBD_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_BD)){
  gene <- geno_BD$gene[i]
  alleles <- geno_BD$GENOTYPED_ALLELES[i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBD_GERM) %in%  paste(gene,alleles,sep='*')
  TRBD_GERM.NEW <- c(TRBD_GERM.NEW,TRBD_GERM[IND])
}

## Remove from TRBJ_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBJ_GERM),'*',fixed=T),'[',1) %in%  geno_BJ$gene)
TRBJ_GERM.NEW <- TRBJ_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_BJ)){
  gene <- geno_BJ$gene[i]
  alleles <- geno_BJ$GENOTYPED_ALLELES[i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBJ_GERM) %in%  paste(gene,alleles,sep='*')
  TRBJ_GERM.NEW <- c(TRBJ_GERM.NEW,TRBJ_GERM[IND])
}


### CHECK IF THE REPLACEMENT IS CORRECT

## Combine the genotyped and others and write to a fasta file for reference
write.fasta(sequences=as.list(gsub(TRBV_GERM.NEW,pattern = '.',replacement = '',fixed = T)),names=names(TRBV_GERM.NEW), v_personal_ref,open="w")
write.fasta(sequences=as.list(gsub(TRBD_GERM.NEW,pattern = '.',replacement = '',fixed = T)),names=names(TRBD_GERM.NEW), d_personal_ref, open="w")
write.fasta(sequences=as.list(gsub(TRBJ_GERM.NEW,pattern = '.',replacement = '',fixed = T)),names=names(TRBJ_GERM.NEW), j_personal_ref, open="w")

# create personal reference db
system(paste(makeblastdb_path, '-parse_seqids -dbtype nucl -in', v_personal_ref,'-out',personal_db_V))
system(paste(makeblastdb_path, '-parse_seqids -dbtype nucl -in', d_personal_ref,'-out',personal_db_D))
system(paste(makeblastdb_path, '-parse_seqids -dbtype nucl -in', j_personal_ref,'-out',personal_db_J))


####################################################################################################################
#################### THIRD IGBLAST + MAKEDB RUNNING ACCORDING TO PERSONAL GENOTYPE REFERENCE #######################
####################################################################################################################

print("RUNNING IGBLAST ON THE PERSONAL REFERENCE")
system(paste(load_ig_environment, igblastn_path, "-germline_db_V", personal_db_V,
             "-germline_db_J", personal_db_J, "-germline_db_D", personal_db_D,
             igblast_constant_args, "-query", novel_igblast_input, "-out", personal_igblast_output))


print("RUNNING MAKEDB ON THE PERSONAL REFERENCE")
write.fasta(sequences = as.list(c(TRBV_GERM.NEW,TRBJ_GERM.NEW)),names = c(names(TRBV_GERM.NEW), names(TRBJ_GERM.NEW)), personal_repo, open="w")
personal_makedb_output <- paste0(sample_path, SAMP, "_genotyped.tab")
system(paste( makedb_com, "-i", personal_igblast_output,
              "-s", novel_igblast_input, '-r', personal_repo, "-o", personal_makedb_output, makedb_additional_parameters))



####################################################################################################################
########################################### WRITE GENOTYPE FILE ####################################################
####################################################################################################################

DATA <- read.delim(personal_makedb_output, header=T, sep = "\t", stringsAsFactors = F)

v_table <- DATA[!grepl(",", DATA$v_call),]
v_table <- v_table %>% dplyr::group_by(v_call) %>% dplyr::summarise(N = n())
v_table$gene <- sapply(strsplit(v_table$v_call, "*", fixed = T), "[", 1)
v_table$allele <- sapply(strsplit(v_table$v_call, "*", fixed = T), "[", 2)
geno_BV$Freq_by_Seq <- unlist(lapply(geno_BV$gene, function(g) {
  if (!g %in% v_table$gene) {
    return("0")
  }
  counts <- v_table$N[v_table$gene==g];
  counts <- sort(counts, decreasing = T);
  return(paste(counts, collapse = ";"))
}))


d_table <- DATA[(!grepl(pattern = ',',DATA$d_call) & DATA$d_call !='None') & (DATA$d_sequence_end - DATA$d_sequence_start >= 8),]
d_table <- d_table[complete.cases(d_table$sequence_id),]
d_table <- d_table %>% dplyr::group_by(d_call) %>% dplyr::summarise(N = n())
d_table$gene <- sapply(strsplit(d_table$d_call, "*", fixed = T), "[", 1)
d_table$allele <- sapply(strsplit(d_table$d_call, "*", fixed = T), "[", 2)
geno_BD$Freq_by_Seq <- unlist(lapply(geno_BD$gene, function(g) {
  if (!g %in% d_table$gene) {
    return("0")
  }
  counts <- d_table$N[d_table$gene==g];
  counts <- sort(counts, decreasing = T);
  return(paste(counts, collapse = ";"))
}))

j_table <- DATA[!grepl(",", DATA$j_call),]
j_table <- j_table %>% dplyr::group_by(j_call) %>% dplyr::summarise(N = n())
j_table$gene <- sapply(strsplit(j_table$j_call, "*", fixed = T), "[", 1)
j_table$allele <- sapply(strsplit(j_table$j_call, "*", fixed = T), "[", 2)
geno_BJ$Freq_by_Seq <- unlist(lapply(geno_BJ$gene, function(g) {
  if (!g %in% j_table$gene) {
    return("0")
  }
  counts <- j_table$N[j_table$gene==g];
  counts <- sort(counts, decreasing = T);
  return(paste(counts, collapse = ";"))
}))


geno <- rbind(geno_BV,geno_BD)
geno <- rbind(geno,geno_BJ)
geno$Freq_by_Clone <- geno$Freq_by_Seq
write.table(geno, file = paste0(genotypes_path, SAMP, "_genotype.tsv"),quote = F,row.names = F,sep = "\t")


####################################################################################################################
############################################# Haplotype Inference ##################################################
####################################################################################################################

hetero_j16 <- "TRBJ1-6" %in% geno$gene
if (hetero_j16) {
  hetero_j16 <- grepl(",", geno$GENOTYPED_ALLELES[geno$gene == "TRBJ1-6"])
}

hetero_d2 <- "TRBD2" %in% geno$gene
if (hetero_d2) {
  hetero_d2 <- grepl(",", geno$GENOTYPED_ALLELES[geno$gene == "TRBD2"])
}

if (!((hetero_j16) | (hetero_d2))) {
  quit()
}

DATA <- read.delim(personal_makedb_output, header=T, sep = "\t", stringsAsFactors = F)
DATA$subject <- SAMP

# Filtering
DATA <- DATA[!grepl(",", DATA$v_call),] # V single assignment

# zero mutations over the V
v_seqs <- sapply(1:nrow(DATA),function(x) substr(DATA$sequence_alignment[x],1,DATA$v_germline_end[x]))
DATA$v_mut <- unlist(tigger::getMutCount(v_seqs, DATA$v_call, germline_db = TRBV_GERM))
DATA <- DATA[DATA$v_mut <= 1,]

del_genes <- geno$gene[grepl("Del", geno$GENOTYPED_ALLELES)]

if (hetero_j16) {
  # Only rearrangments with single assignment of TRBJ1-6
  DATA_J <- DATA[grepl("J1-6", DATA$j_call),]
  DATA_J <- DATA_J[!grepl(",", DATA_J$j_call),]
  DATA_J$J_SEQ_LENGTH <- DATA_J$j_sequence_end - DATA_J$j_sequence_start + 1
  DATA_J <- DATA_J[DATA_J$J_SEQ_LENGTH > 10,]
  
  TRBJ1_6_01_REA <- nrow(DATA_J[DATA_J$j_call == "TRBJ1-6*01",])
  total <- nrow(DATA_J)
  
  if ((TRBJ1_6_01_REA / total > 0.3) & (TRBJ1_6_01_REA / total < 0.7)) {
    haplo_j1_6 <- createFullHaplotype(DATA_J,toHap_col = "v_call", hapBy_col = "j_call", chain = "TRB", hapBy = "TRBJ1-6", toHap_GERM = TRBV_GERM, rmPseudo = F)
    haplo_j1_6$TOTAL <- total
    
    if (length(del_genes) > 0) {
      for (gene in del_genes) {
        haplo_j1_6[haplo_j1_6$gene == gene, c(3:5, 9)] <- c("Del", "Del", "Del", 1000)
      }
    }
    
    write.table(haplo_j1_6, file = paste0(sample_path, SAMP, "_haplo_J1_6.tab"), quote = F, row.names = F, sep = "\t")
  }
}


if (hetero_d2) {
  # Only rearrangments with single assignment of TRBD2
  DATA$d_germline_length <- as.integer(DATA$d_germline_end) - as.integer(DATA$d_germline_start) + 1
  DATA_D_geno <- DATA[(!grepl(pattern = ',',DATA$d_call) & DATA$d_call!='None') & (DATA$d_germline_length >= 9),]
  DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno$sequence_id),]
  DATA_D_geno <- DATA_D_geno[grepl("D2", DATA_D_geno$d_call),]
  
  # filter by zero mutations over the D segment
  # extract d sequence in the direct orientation
  DATA_D_reg <- DATA_D_geno[DATA_D_geno$d_germline_start < DATA_D_geno$d_germline_end,]
  DATA_D_reg$d_seq <- substr(DATA_D_reg$sequence, DATA_D_reg$d_sequence_start, DATA_D_reg$d_sequence_end)
  
  # extract convert d sequence in the inverted orientation to the direct orientation
  DATA_D_inv <- DATA_D_geno[DATA_D_geno$d_germline_start > DATA_D_geno$d_germline_end,]
  DATA_D_inv$d_seq <- substr(DATA_D_inv$sequence, DATA_D_inv$d_sequence_start, DATA_D_inv$d_sequence_end)
  DATA_D_inv$d_seq <- stringi::stri_reverse(DATA_D_inv$d_seq)
  DATA_D_inv$d_seq <- gsub("A", "t", DATA_D_inv$d_seq)
  DATA_D_inv$d_seq <- gsub("T", "a", DATA_D_inv$d_seq)
  DATA_D_inv$d_seq <- gsub("G", "c", DATA_D_inv$d_seq)
  DATA_D_inv$d_seq <- gsub("C", "g", DATA_D_inv$d_seq)
  DATA_D_inv$d_seq <- toupper(DATA_D_inv$d_seq)
  
  d_germ_end <- DATA_D_inv$d_germline_start
  DATA_D_inv$d_germline_start <- DATA_D_inv$d_germline_end
  DATA_D_inv$d_germline_end <- d_germ_end
  
  DATA_D_geno <- rbind(DATA_D_reg, DATA_D_inv)
  
  DATA_D_geno$mut_d <- unlist(lapply(1:nrow(DATA_D_geno), function(i) {
    mut <- 0;
    row_seq <- unlist(strsplit(DATA_D_geno$d_seq[[i]], ""));
    allele_seq <- unlist(strsplit(TRBD_GERM[[DATA_D_geno$d_call[[i]]]], ""));
    for (pos in DATA_D_geno$d_germline_start[[i]]:(DATA_D_geno$d_germline_end[[i]])) {
      if (row_seq[pos - (DATA_D_geno$d_germline_start[[i]] - 1)] != allele_seq[pos]) {
        mut <- mut + 1;
      }
    }
    mut
  }))
  
  DATA_D_geno <- DATA_D_geno[DATA_D_geno$mut_d == 0,]
  TRBD2_01_REA <- nrow(DATA_D_geno[DATA_D_geno$d_call == "TRBD2*01",])
  total <- nrow(DATA_D_geno)
  
  if ((TRBD2_01_REA / total > 0.3) & (TRBD2_01_REA / total < 0.7)) {
    haplo_d2 <- createFullHaplotype(DATA_D_geno, toHap_col = "v_call", hapBy_col = "d_call", chain = "TRB", hapBy = "TRBD2", toHap_GERM = TRBV_GERM, rmPseudo = F, kThreshDel = 5)
    haplo_d2$TOTAL <- total
    
    if (length(del_genes) > 0) {
      for (gene in del_genes) {
        haplo_d2[haplo_d2$gene == gene, c(3:5, 9)] <- c("Del", "Del", "Del", 1000)
      }
    }
    
    write.table(haplo_d2, file = paste0(sample_path, SAMP, "_haplo_D2.tab"), quote = F, row.names = F, sep = "\t")
  }
}
