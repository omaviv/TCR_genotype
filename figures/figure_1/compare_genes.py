# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 10:52:52 2019

@author: aviv
"""
import sys

import multiprocessing
from multiprocessing import Pool
from functools import partial

import numpy as np
import pandas as pd
from .patterns_by_primers import *


def compoare_imgt_sequences(seq1, seq2):
	distance = 0
	for i in range(min(len(seq1), len(seq2))):
		if not ((seq1[i] == 'n') | (seq2[i] == 'n') | (seq1[i] == seq2[i])):
			distance += 1
	
	return distance


def compare_gene_allele_sequences(gene_alleles):
	genes = list(gene_alleles.keys())
	distances = {}
	
	for gene in genes:
		distances[gene] = {gene:0}
		
	for i in range(len(genes)-1):
		gene1 = genes[i]
		
		for j in range(i+1, len(genes)):
			distance = 500			
			gene2 = genes[j]
			
			for allele_seq1 in gene_alleles[gene1].values():
				for allele_seq2 in gene_alleles[gene2].values():
					distance = min(distance, compoare_imgt_sequences(allele_seq1, allele_seq2))
				
			distances[gene1][gene2] = distance
			distances[gene2][gene1] = distance
	
	return distances


# check the global score according to the lower length between the allele known sequence lengths
# (not any more) and the length of 322 nt.
def global_alignment(allele_seq1, allele_seq2, reverse = False):
	seq_len = min(len(allele_seq1), len(allele_seq2)) + 1
	seq1 = allele_seq1[:seq_len]
	seq2 = allele_seq2[:seq_len]

	seq1 = seq1.replace(".", "")	
	seq2 = seq2.replace(".", "")
	seq_len1 = len(seq1) + 1
	seq_len2 = len(seq2) + 1
	m = np.zeros([seq_len1 , seq_len2])

	indel_penlty = 1	
	mismatch_penlty = 1	
	
	if reverse:
		seq1 = seq1[::-1]
		seq2 = seq2[::-1]
    
	for i in range(1, seq_len1):
		for j in range(1, seq_len2):
			insertion = m[i-1, j] + indel_penlty
			deletion = m[i, j-1] + indel_penlty
			match = m[i-1, j-1]
			if (seq1[i-1] != seq2[j-1]) & (seq1[i-1] != 'n') & ('n' != seq2[j-1]):
				match += mismatch_penlty
			
			m[i, j] = min(insertion, deletion, match)
	
	return m[-1, -1]
	

def compare2genes(args, reverse, gene_alleles):
	gene1 = args["gene1"] 
	gene2 = args["gene2"]
	
	
	distance = 500
	for allele_seq1 in gene_alleles[gene1].values():
		for allele_seq2 in gene_alleles[gene2].values():
			distance = min(distance, global_alignment(allele_seq1, allele_seq2, reverse))
			
	return distance


def compare_gene_sequences(gene_alleles, reverse = False):
	genes = list(gene_alleles.keys())
	distances = {}	

	for gene in genes:
		distances[gene] = {gene:0}
		
	gene2compare = [{"gene1":genes[i],"gene2":genes[j]}
	 for i in range(len(genes)-1) for j in range(i+1, len(genes))]

	with Pool(6) as p:
		fun_and_arg = partial(compare2genes, reverse=reverse, gene_alleles=gene_alleles)
		distance_values = p.map(fun_and_arg, gene2compare)
	
	for i in range(len(gene2compare)):
		gene1 = gene2compare[i]["gene1"]
		gene2 = gene2compare[i]["gene2"]
		distance = distance_values[i]
		
		distances[gene1][gene2] = distance
		distances[gene2][gene1] = distance
	
	return distances


def extract_gene_alleles(allele_dict):
	gene_names = allele_dict.keys()
	gene_names = [allele.split("*")[0] for allele in gene_names]
	gene_names = set(gene_names)
	gene_names = list(gene_names)
		
	temp = list(sorted(list(allele_dict.keys()))[0].split("*")[1])
	allele_naming = "*"
	for c in temp:
		allele_naming += c
		if c == "0":
			break		
	
	gene_seqs = {}	
	for gene in gene_names:
		gene_seqs[gene] = {}
		
		for i in range(1, 10):
			allele_name = gene + allele_naming + str(i)
			if allele_name not in alleles.keys():
				break
			
			gene_seqs[gene][allele_name] = alleles[allele_name]							
	
	return gene_seqs



if __name__ == "__main__":
	project_path = "tcr_genotype/"
	required_files_path = project_path + "figures/data/"
	figure_path = project_path + "figures/figure_1/"
	
	# Gene distances according to the full V sequences
	# alleles_file = required_files_path + "TRBV_all_genes.fasta"
	# output_file = required_files_path + "TRBV_GENE_DISTANCES.tab"
	
	# Gene distances according to BIOMED-2 coverage
	alleles_file = required_files_path + "BIOMED2_TRBV_all_genes.fasta"
	output_file = required_files_path + "BIOMED2_TRBV_GENE_DISTANCES.tab"
	
	# Gene distances according to Adaptive Biotechnologies coverage
	# alleles_file = required_files_path + "Adaptive_TRBV_all_genes.fasta"
	# output_file = required_files_path + "Adaptive_TRBV_GENE_DISTANCES.tab"
	
	alleles = read_fasta(alleles_file)

	# For alignment by the shared gene sequence
	genes = extract_gene_alleles(alleles)
 					
	distances = compare_gene_sequences(genes)
	distances_df = pd.DataFrame.from_dict(distances)
	distances_df['GENES'] = distances_df.index
	cols = distances_df.columns.tolist()
	cols = cols[-1:] + cols[:-1]
	distances_df = distances_df[cols]
	distances_df.to_csv(output_file, sep = "\t", index=False)	
 		
