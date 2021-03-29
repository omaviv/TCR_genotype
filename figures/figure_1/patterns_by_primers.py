# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 14:13:14 2019

@author: aviv
"""
import re, xlrd
import pandas as pd

def get_primers_by_xls(filename):
	book = xlrd.open_workbook(filename)
	sheet = book.sheet_by_index(0)
	primers = []
	
	for row in range(1, sheet .nrows):
	    primers.append(sheet.cell(row, 1).value.strip().lower())
	
	primers = set(primers)
	primers = list(primers)
	
	return primers


def read_fasta(filename):
	head_pattern = re.template("^>")
	name = None
	seq = ""	
	alleles = {}
	fastafile = open(filename, "r")	
	
	for line in fastafile:
		if re.search(head_pattern, line):
			if name:
				alleles[name]=seq.lower()
			name = line.replace(">","").strip().replace("D", "")
			seq = ""
		else:
			seq += line.strip()
	
	alleles[name]=seq.lower()
	fastafile.close()
	
	return alleles


def extract_genes(allele_dict):
	gene_names = allele_dict.keys()
	gene_names = [allele.split("*")[0] for allele in gene_names]
	gene_names = set(gene_names)
	gene_names = list(gene_names)
	
	gene_seqs = {}
	
	temp = list(sorted(list(allele_dict.keys()))[0].split("*")[1])
	allele_naming = "*"
	for c in temp:
		allele_naming += c
		if c == "0":
			break		
	
	all_alleles = sorted(list(alleles.keys()))
	gene = ""
	for allele in all_alleles:
		if allele == all_alleles[0]:
			gene = allele.split("*")[0]
			compared_allele = allele
			seq = list(allele_dict[compared_allele])
			continue
		
		elif gene != allele.split("*")[0]:
			gene_seqs[gene] = "".join(seq)
			gene = allele.split("*")[0]
			compared_allele = allele
			seq = list(allele_dict[compared_allele])
			continue
		
		temp_seq = list(allele_dict[allele])
		for i in range(min(len(seq), len(temp_seq))):
				if (seq[i] != temp_seq[i]) & (temp_seq[i] != "."):
					seq[i] = "n"
			
	gene_seqs[gene] = "".join(seq)
	
	return gene_seqs


def primers_pos(primers, genes):
	"""
	This function gets dictionary of 
	"""
	primer_pos = {}
	primer_ind = 1
	for primer in primers:
		primer_pos[primer]={}
		print("primer " + str(primer_ind) + " out of " + str(len(primers)) + " primers")
		primer_ind+=1
		for s in range(0, len(primer)-3):
			for i in range(s+1, len(primer)-2):
				for j in range(i+1, len(primer)-1):
					for k in range(j+1, len(primer)):
						primer_pattern = list(primer)
						primer_pattern[s] = 'n'
						primer_pattern[i] = 'n'
						primer_pattern[j] = 'n'
						primer_pattern[k] = 'n'
						primer_pattern = "".join(primer_pattern)					
						
						primer_pattern = primer_pattern.replace("a","[a]")
						primer_pattern = primer_pattern.replace("t","[t]")
						primer_pattern = primer_pattern.replace("g","[g]")
						primer_pattern = primer_pattern.replace("c","[c]")
						primer_pattern = primer_pattern.replace("n","[a-z]")
						
						primer_pattern = primer_pattern.replace("][","][.]*[")
						
						for gene in genes.keys():
							temp = re.search(primer_pattern, genes[gene])
							if re.search(primer_pattern, genes[gene]):
								if gene in primer_pos[primer].keys():
									primer_pos[primer][gene] = max(temp.span()[1],primer_pos[primer][gene])
								else:
									primer_pos[primer][gene] = temp.span()[1]
			
	return primer_pos


def collapse_patterns(genes, alleles, pat = "ap"):
	new_patterns = {}
	patterns_alleles = {}
	for gene in genes:
		taken = []
		pattern_index = 0
		for i in range(1, 10):
			if not i in taken:
				compared_allele = gene + "*0" + str(i)
				if compared_allele not in alleles.keys():
					break
				
				pattern_index += 1
				taken.append(i)
				pattern_name = gene + "*" + pat + "0" + str(pattern_index)
				patterns_alleles[pattern_name] = ["0" + str(i)]				
				new_seq = alleles[compared_allele]
				
				for j in range(i+1, 10):
					if not j in taken:
						temp_allele = gene + "*0" + str(j)
						if temp_allele not in alleles.keys():
							break
										
						if alleles[temp_allele] in new_seq:
							taken.append(j)
							patterns_alleles[pattern_name].append("0" + str(j))
						elif new_seq in alleles[temp_allele]:
							taken.append(j)
							patterns_alleles[pattern_name].append("0" + str(j))
							new_seq = alleles[temp_allele]
				
				new_patterns[pattern_name] = new_seq
	
	return new_patterns, patterns_alleles
				

def full_pattern_sequence(seq1, seq2):
	if len(seq1) < len(seq2):
		temp = seq1
		seq1 = seq2
		seq2 = temp
	
	for i in range(len(seq2)):
		if (seq1[i] != seq2[i]) & (seq2[i] != "."):
			seq1[i] = 'n'



if __name__ == "__main__":
	project_path = "tcr_genotype/"
	required_files_path = project_path + "figures/data/"
	figure_path = project_path + "figures/figure_1/"
	
	primers_file = figure_path + "BIOMED2_V_PRIMERS.xls"
	reference_file = required_files_path + "TRBV_all_genes.fasta" # for example IGHV.TRBV
	output_reference_file = "BIOMED2_TRBV_all_genes.fasta"
	output_table = figure_path + "BIOMED2_TRBV.tab"
	family_mismatches_output = figure_path + "BIOMED2_TRBV_mismatches.tab"
	allele_prefix = "bp" # the prefix before the allele number
	
	
	primers = get_primers_by_xls(primers_file)
	alleles = read_fasta(reference_file)
	genes = extract_genes(alleles)
	
	# find the primer last positiion over the genes sequences
	gene_p_pos = primers_pos(primers, genes)
	
	# find the position to cut each gene.
	gene_pos2cut = {} 
	for gene in genes.keys():
		gene_pos2cut[gene] = 0
		for primer in primers:
			if gene in gene_p_pos[primer]:
				gene_pos2cut[gene] = max(gene_pos2cut[gene], gene_p_pos[primer][gene])
		
		if gene_pos2cut[gene] == 0:
			print(gene)
	
	
	biomed_alleles = {}
	for allele in alleles.keys():
		if gene_pos2cut[allele.split("*")[0]]:
			gene = allele.split("*")[0]
			gene_seq2primer = "." * (gene_pos2cut[gene] + 1)
			biomed_alleles[allele] = gene_seq2primer + alleles[allele][gene_pos2cut[gene]+1:]
	
	biomed_alleles, patterns_alleles = collapse_patterns(genes=genes, alleles=biomed_alleles, pat = allele_prefix)
	
	# create new alleles file
	alleles_names = list(biomed_alleles.keys())
	alleles_names = sorted(alleles_names)
	new_ref_file = open(output_reference_file,"w")
	for allele in alleles_names:
		new_ref_file.write(">" + allele + "\n")
		new_ref_file.write(biomed_alleles[allele] + "\n")
	new_ref_file.close()
	
	names = list(patterns_alleles.keys())
	names = sorted(names)
	for pat in names:
		patterns_alleles[pat] = ", ".join(patterns_alleles[pat])
		print(pat + ": " + patterns_alleles[pat])
	
	
	patterns_df = pd.DataFrame.from_dict(patterns_alleles, orient = 'index')
	patterns_df['PATTERN'] = patterns_df.index
	patterns_df.reset_index()
	patterns_df[['GENE', 'PATTERN']] = patterns_df['PATTERN'].str.split("*", expand=True)
	patterns_df['ALLELES'] = patterns_df[0]
	patterns_df = patterns_df[['GENE', 'PATTERN', 'ALLELES']]
	patterns_df = patterns_df.sort_values(['GENE', 'PATTERN'])
	patterns_df.to_csv(output_table, sep = "\t", index=False)
	
	
	family_mismatches = {}
	families = set([allele_name.split("*")[0].split("-")[0] for allele_name in biomed_alleles.keys()])
	for family in families:
		family_alleles = [allele_name for allele_name in biomed_alleles.keys() if family == allele_name.split("*")[0].split("-")[0]]
		
		if len(family_alleles) == 1:
			continue
		
		for allele_1 in family_alleles:
			allele_seq_1 = biomed_alleles[allele_1]			
			for allele_2 in family_alleles:
				if allele_2 == allele_1:
					continue
				
				key = allele_1 + "," + allele_2
				
				family_mismatches[key] = []
				allele_seq_2 = biomed_alleles[allele_2]			
				for i in range(min(len(allele_seq_1), len(allele_seq_2))):
					nuc_1 = allele_seq_1[i]
					nuc_2 = allele_seq_2[i]
					
					if (nuc_1 != nuc_2) & (nuc_1 != ".") & (nuc_2 != "."):
						family_mismatches[key].append(nuc_1.upper() + str(i+1) + nuc_2.upper())
					
				
				family_mismatches[key] = ",".join(family_mismatches[key])
			
	
	family_mismatches_df = pd.DataFrame.from_dict(family_mismatches, orient = 'index')
	family_mismatches_df.reset_index(level=0, inplace=True)
	family_mismatches_df.rename(columns={"index":"ALLELES", 0:"MISMATCHES"}, inplace=True)
	
	allele_cols = family_mismatches_df["ALLELES"].str.split(',', expand=True)
	family_mismatches_df["ALLELE_1"] = allele_cols[0]
	family_mismatches_df["ALLELE_2"] = allele_cols[1]
	family_mismatches_df.drop(columns=["ALLELES"], inplace = True)
	
	family_mismatches_df = family_mismatches_df[["ALLELE_1", "ALLELE_2", "MISMATCHES"]]
	
	family_mismatches_df.to_csv(family_mismatches_output, sep = "\t", index=False)
	
	