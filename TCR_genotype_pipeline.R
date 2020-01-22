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


####################################################################################################################
######################### PARAMETERS TO INTIAL BEFORE FIRST RUN ####################################################
####################################################################################################################

# The number of threads to use for IgBlast.
num_of_threads <- 44

# IgBlast path. In case that igblast can run simply by igblastn command it could be just "igblastn"
igblastn_path <- "/<path_to_igblast>/igblastn"


# the paths to the germline databases which done by makeblastdb
germline_db_V <- "/<path_to_germline_databases>/<TRBV_database_name>"
germline_db_D <- "/<path_to_germline_databases>/<TRBD_database_name>"
germline_db_J <- "/<path_to_germline_databases>/<TRBJ_database_name>"

# in
load_ig_environment <- "export IGDATA=/home/bcrlab/eitan/data/igblast/; nice -19"

# the path for makeblastdb tool
makeblastdb_path <- "/<Path_to_makeblastdb>/makeblastdb"

# the paths of the fasta files which conatins the germline reference in IMGT format.
TRBV_imgt_ref <- "/<path_to the_ger>/<TRBV>.fasta"
TRBD_imgt_ref <- "/<path_to the_ger>/<TRBD>.fasta"
TRBJ_imgt_ref <- "/<path_to the_ger>/<TRBJ>.fasta"


####################################################################################################################
############################### Handling the flags  ################################################################
####################################################################################################################

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="fasta file name", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="reference file name", metavar="character"),
  make_option(c("-p", "--path"), type="character", default="/localdata/aviv/genotypes", 
              help="output path folder", metavar="character"))

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


if (!(tolower(opt$sequence_type) %in% c("full", "bio", "adapt"))){
  print_help(opt_parser)
  stop("sequence_type must be (full\ bio \ adapt)", call.=FALSE)
}


igblast_input <- opt$file
SAMP <- opt$sample
output_path <- opt$path
sequencing_length <- tolower(opt$sequence_type)

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
tcrb_repo <- paste0(output_path, "/TCRB_repo.fasta") # The reference sequences into the "REPO" file for MakeDB.


### PERSONAL NOVEL IGBLAST
personal_novel_fasta <- paste0(sample_path, "with_novel_V_personal_ref.fasta")
personal_novel_germline_db_V <- paste0(sample_path, "/personal_ref/novel_V")

### PERSONAL GENOTYPE PARAMETERS
personal_ref_path <- paste0(sample_path, "personal_ref/")
v_personal_ref <- paste0(personal_ref_path, "V_personal_ref.fasta")
d_personal_ref <- paste0(personal_ref_path, "D_personal_ref.fasta")
j_personal_ref <- paste0(personal_ref_path, "J_personal_ref.fasta")
persona_repo <- paste0(sample_path, "personal_repo.fasta")
personal_fasta <- paste0(sample_path, SAMP, ".fasta")

personal_db_V <- paste0(personal_ref_path, "genotyped_BV")
personal_db_D <- paste0(personal_ref_path, "genotyped_BD")
personal_db_J <- paste0(personal_ref_path, "genotyped_BJ")

personal_igblast_output <- paste0(sample_path, SAMP, "_genotyped.fmt7")

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
names(TRBV_GERM) <- sapply(strsplit(names(TRBV_GERM),"|",fixed = T),"[",2) 
TRBV_GERM <- toupper(unlist(TRBV_GERM))

TRBD_GERM <- read.fasta(TRBD_imgt_ref, as.string = TRUE)
names(TRBD_GERM) <- sapply(strsplit(names(TRBD_GERM),"|",fixed = T),"[",2) 
TRBD_GERM <- toupper(unlist(TRBD_GERM))

TRBJ_GERM <- read.fasta(TRBJ_imgt_ref, as.string = TRUE)
names(TRBJ_GERM) <- sapply(strsplit(names(TRBJ_GERM),"|",fixed = T),"[",2) 
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

# igblast command
system(paste(load_ig_environment, igblastn_path,"-germline_db_V", germline_db_V, 
             "-germline_db_D", germline_db_D, "-germline_db_J", germline_db_J, igblast_constant_args,
             "-query", igblast_input, "-out", igblast_output))

print("RUNNING MAKEDB ON THE FIRST IGBLAST OUTPUT")
makedb_output <- paste0(makedb_out_dir, SAMP, ".tab")
system(paste(load_makedb_env, makedb_com, "-i", igblast_output, "-s", igblast_input, '-r', tcrb_repo, "-o", makedb_output))


####################################################################################################################
####################################### SEARCHING FOR NOVEL ALLELLES ###############################################
####################################################################################################################
# read the table after the first makedb
DATA <- read.delim(makedb_output, header=T, sep = "\t", stringsAsFactors = F)

print("SEARCHING FOR NOVEL ALLELLES")
novel_df_H <- findNovelAlleles(DATA, TRBV_GERM, germline_min = 1, min_seqs = 1, min_frac = 0.35)
#Remove ambigues assignments for the same sequences for two differnet novel allels
novel_df_H <- novel_df_H %>% group_by(NOVEL_IMGT) %>% slice(which.max(PERFECT_MATCH_COUNT/GERMLINE_CALL_COUNT))

#Remove novel allelles that are of the same distance from two existent alleles (leave)
novel_df_H$MUTDIST <- apply(novel_df_H,1,function(x){sum(s2c(as.character(x['NOVEL_IMGT']))!=s2c(as.character(x['GERMLINE_IMGT'])))})
novel_df_H <- novel_df_H %>%  filter( !is.na(POLYMORPHISM_CALL)) %>% mutate(GENE= sapply(strsplit(GERMLINE_CALL,'*',fixed=T),'[',1))
if(nrow(novel_df_H)!=0) {
  new_novel_df_H <- novel_df_H %>%  group_by(GENE,MUTDIST) %>% slice(1) %>% ungroup
} else {
  new_novel_df_H <- novel_df_H
}

## ADD novel alleles to the genotype database
novel_genotype <- new_novel_df_H$NOVEL_IMGT
novel_genotype <- sapply(novel_genotype,function(x){gsub('-','.',x,fixed = T)})
names(novel_genotype) <- new_novel_df_H$POLYMORPHISM_CALL
for (novel_allele in names(novel_genotype)) {
  if (!(novel_allele %in% names(TRBV_GERM))) {
    TRBV_GERM <- c(TRBV_GERM,novel_genotype[novel_allele])
  }
}

## Combine the genotyped and others and write to a fasta file for reference
write.fasta(sequences=as.list(gsub(TRBV_GERM, pattern = '.', replacement = '', fixed = T)),
            names=names(TRBV_GERM), personal_novel_fasta, open="w")

## Write the VDJ sequences to realign
info.strt.idx <- which(names(DATA) == 'JUNCTION')
seq.names <- sapply(1:nrow(DATA),function(x){ paste0(names(DATA)[c(1,(info.strt.idx+1):(ncol(DATA)))],
                                                     rep('=',length(c(1,(info.strt.idx+1):(ncol(DATA))))),
                                                     DATA[x,c(1,(info.strt.idx+1):(ncol(DATA)))],
                                                     collapse = '|')})
seq.names <- gsub('SEQUENCE_ID=', '', seq.names, fixed = T)
write.fasta(sequences=as.list(DATA$SEQUENCE_VDJ), names=seq.names, paste0(SAMP, ".fasta"), open="w")

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
igblast_input <- paste0(pathname, "/", SAMP, ".fasta")
igblast_output <- paste0(pathname,"/igblast_novel.out")
system(paste(load_ig_environment, igblastn_path, "-germline_db_V ", personal_novel_germline_db_V, 
             "-germline_db_D", germline_db_D, "-germline_db_J", germline_db_J, 
             igblast_constant_args, "-query", igblast_input, "-out", igblast_output))


print("RUNNING MAKEDB ON THE IGBLAST OUTPUT")
novel_repo <- paste0(output_path, "/", SAMP, "/", "novel_repo.fasta")
write.fasta(sequences = as.list(c(TRBV_GERM, TRBJ_GERM)), names = c(names(TRBV_GERM), names(TRBJ_GERM)),
            novel_repo, open="w")


print("RUNNING MAKEDB ON THE 2ND IGBLAST OUTPUT")

makedb_output <- paste0(sample_path, SAMP, "_novel.tab")
system(paste(load_makedb_env, makedb_com, "-i",  igblast_output, "-s", igblast_input, "-r", novel_repo, "-o", makedb_output))

####################################################################################################################
########################################## GENOTYPE INFERRING ######################################################
####################################################################################################################
DATA <- read.delim(makedb_output, header=T, sep = "\t", stringsAsFactors = F)

DATA_V_SA <- DATA[!grepl(pattern = ',',DATA$V_CALL),]
geno_H <- inferGenotypeBayesian(DATA_V_SA, germline_db = TRBV_GERM, find_unmutated = TRUE,
                                novel = new_novel_df_H,v_call = 'V_CALL')
geno_H$GENOTYPED_ALLELES <-apply(geno_H[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:5]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})

DATA_D_geno <- DATA[(!grepl(pattern = ',',DATA$D_CALL) & DATA$D_CALL!='None') & (DATA$D_GERM_LENGTH >= 9),]
DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno$SEQUENCE_ID),]

genoD_H <- inferGenotypeBayesian(DATA_D_geno, find_unmutated = FALSE, germline_db = TRBD_GERM, v_call = 'D_CALL')
genoD_H$GENOTYPED_ALLELES <-apply(genoD_H[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})

DATA_J_SA <- DATA[!grepl(pattern = ',',DATA$J_CALL),]
genoJ_H <- inferGenotypeBayesian(DATA, germline_db = TRBJ_GERM, find_unmutated = FALSE, v_call = 'J_CALL')
genoJ_H$GENOTYPED_ALLELES <-apply(genoJ_H[,c(2,6:9)],1,function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]),','))[1:m],collapse=",")})


geno <- rbind(geno_H,genoD_H)
geno <- rbind(geno,genoJ_H)
write.table(geno, file = paste0(genotypes_path, SAMP, "_geno_H.tab"),quote = F,row.names = F,sep = "\t")


####################################################################################################################
#################################### CREATE PERSONAL GENOTYPE REFERENCE ############################################
####################################################################################################################

## Remove from TRBV_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBV_GERM),'*',fixed=T),'[',1) %in%  geno_H$GENE)
TRBV_GERM.NEW <- TRBV_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_H)){
  
  gene <- geno_H$GENE[i]
  
  alleles <- geno_H$GENOTYPED_ALLELES[i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBV_GERM) %in%  paste(gene,alleles,sep='*')
  TRBV_GERM.NEW <- c(TRBV_GERM.NEW,TRBV_GERM[IND])
}


## Remove from TRBD_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBD_GERM),'*',fixed=T),'[',1) %in%  genoD_H$GENE)
TRBD_GERM.NEW <- TRBD_GERM[NOTGENO.IND]

for(i in 1:nrow(genoD_H)){
  gene <- genoD_H$GENE[i]
  alleles <- genoD_H$GENOTYPED_ALLELES[i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBD_GERM) %in%  paste(gene,alleles,sep='*')
  TRBD_GERM.NEW <- c(TRBD_GERM.NEW,TRBD_GERM[IND])
}

## Remove from TRBJ_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBJ_GERM),'*',fixed=T),'[',1) %in%  genoJ_H$GENE)
TRBJ_GERM.NEW <- TRBJ_GERM[NOTGENO.IND]

for(i in 1:nrow(genoJ_H)){
  gene <- genoJ_H$GENE[i]
  alleles <- genoJ_H$GENOTYPED_ALLELES[i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBJ_GERM) %in%  paste(gene,alleles,sep='*')
  TRBJ_GERM.NEW <- c(TRBJ_GERM.NEW,TRBJ_GERM[IND])
}


### CHECK IF THE REPLACEMENT IS CORRECT

## Combine the genotyped and others and write to a fasta file for reference
write.fasta(sequences=as.list(gsub(TRBV_GERM.NEW,pattern = '.',replacement = '',fixed = T)),names=names(TRBV_GERM.NEW), v_personal_ref,open="w")
write.fasta(sequences=as.list(gsub(TRBD_GERM.NEW,pattern = '.',replacement = '',fixed = T)),names=names(TRBD_GERM.NEW), d_personal_ref, open="w")
write.fasta(sequences=as.list(gsub(TRBJ_GERM.NEW,pattern = '.',replacement = '',fixed = T)),names=names(TRBJ_GERM.NEW), j_personal_ref, open="w")


## Write the VDJ sequences to realign
info.strt.idx <- which(names(DATA) == 'JUNCTION')
seq.names <- sapply(1:nrow(DATA),function(x){ paste0(names(DATA)[c(1,(info.strt.idx+1):(ncol(DATA)-1))],
                                                     rep('=',length(c(1,(info.strt.idx+1):(ncol(DATA)-1)))),
                                                     DATA[x,c(1,(info.strt.idx+1):(ncol(DATA)-1))],
                                                     collapse = '|')})
seq.names <- gsub('SEQUENCE_ID=','',seq.names,fixed = T)

## Write the VDJ sequences to realign
write.fasta(sequences=as.list(DATA$SEQUENCE_VDJ),names=seq.names,personal_fasta,open="w")

print("CREATING A PERSONAL REFERENCE DATABASE")
system(paste(makeblastdb_path, '-parse_seqids -dbtype nucl -in', v_personal_ref,'-out',personal_db_V))
system(paste(makeblastdb_path, '-parse_seqids -dbtype nucl -in', d_personal_ref,'-out',personal_db_D))
system(paste(makeblastdb_path, '-parse_seqids -dbtype nucl -in', j_personal_ref,'-out',personal_db_J))


####################################################################################################################
#################### THIRD IGBLAST + MAKEDB RUNNING ACCORDING TO PERSONAL GENOTYPE REFERENCE #######################
####################################################################################################################

print("RUNNING IGBLAST ON THE PERSONAL REFERENCE")
system(paste(load_ig_environment, igblastn_path, "-germline_db_V", personal_db_V, 
             "-germline_db_J ",personal_db_J, "-germline_db_D", personal_db_D, 
             igblast_constant_args, "-query", igblast_input, "-out", personal_igblast_output))


print("RUNNING MAKEDB ON THE PERSONAL REFERENCE")
write.fasta(sequences = as.list(c(TRBV_GERM.NEW,TRBJ_GERM.NEW)),names = c(names(TRBV_GERM.NEW), names(TRBJ_GERM.NEW)), persona_repo, open="w")
personal_makedb_output <- paste0(sample_path, SAMP, "_genotyped.tab")
system(paste(load_makedb_env, makedb_com, "-i", personal_igblast_output, 
             "-s", igblast_input, '-r', persona_repo, "-o", personal_makedb_output))

