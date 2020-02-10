# TCR_genotype

## Requirements:
1. IgBlast tool
2. makeblastdb tool
3. MakeDb tool
4. R and the following R packages: seqinr, tigger, optparse, gtools, ggplot2, stringr, dplyr, data.table, alakazam, and reshape2
5. The fasta file refernces of the 3 segments: TRBV, TRBD, and TRBJ in IMGT format (gapped sequences). (the references can be found in the "Script_references" folder)
6. Blast databases for germline TRBV, TRBD, and TRBJ sequences as describde over IgBlast set up: https://ncbi.github.io/igblast/cook/How-to-set-up.html (the references can be found in the "makeblastdb_references" folder)

## Initialize
Before running the pipeline one has to initialize the parameters at the beginning of the code "PARAMETERS TO INTIAL BEFORE FIRST RUN" (The first section after the libraries).
The parameter are described at the same section.

## How to run the script
The command for runnig the script with the essential flags is:
nice -19 Rscript <ptah_to_the_script>/TCR_genotype_pipeline_test.R -f <fasta_file> -s <Sample_name> -p <output_dir>
"output_dir" - in which directory to store all the data

Another optional flag is: "-c <CONSCOUNT>" - Only sequences with atleast the CONSCOUNT value would include to infer the genotype. The default is 1.  

