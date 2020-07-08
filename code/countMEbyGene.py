#!/usr/bin/env python
#
# 
# ****Function Description***
# Counts number of mutation events per gene
#
# Author(s) 
#	D. J. Edwards (David.Edwards@monash.edu) 
#
# Example command:
'''
python countMEbyGene.py -h
'''
#
# Last modified - 09/07/2020
# Change - simplified output (no reporting of codons and sites beyond count)

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nCounts Mutation Events per gene')
	parser.add_argument('-i', '--input_file', type=str, required=True, help='ME file')
	parser.add_argument('-p', '--output_prefix', type=str, required=True, help='Prefix to add to output file')
	return parser.parse_args()

def readInput(input_file_name):
	input_file_handle = open(input_file_name, 'r')
	input_file = input_file_handle.readlines()
	input_file_handle.close()
	return input_file

def writeOutput(output_file_name,output):
	output_file_handle = open(output_file_name,"w")
	output_file_handle.write(output)
	output_file_handle.close()
	return

def countEvents(mutation_events_file,output):
	genes = []
	gene_counts = []
	gene_counts_Syn = []
	gene_counts_NonSyn = []
	gene_counts_Ambiguous = []
	gene_hits = []
	gene_hits_counts = []
	gene_hits_NonSyn_counts = []
	gene_hits_Syn_counts = []
	gene_hits_Amb_counts = []
	gene_Codons =[]
	gene_Codons_counts =[]
	gene_Codons_NonSyn_counts = []
	gene_Codons_Syn_counts = []
	gene_Codons_Amb_counts = []

	for line in mutation_events_file[1:]:
		entry = line.split('\t')
		if entry[6] != '-': #test for intragenic SNP
			gene = entry[6]
			snp = entry[0]
			codon = entry[8]
			if gene in genes:
				index = genes.index(gene)
				gene_counts[index] += 1
				change_type = entry[14]
				if change_type == 'NS':
					gene_counts_NonSyn[index] += 1
					if snp in gene_hits[index]:
						snp_index = gene_hits[index].index(snp)
						gene_hits_counts[index][snp_index] += 1
						gene_hits_NonSyn_counts[index][snp_index] += 1
					else:
						gene_hits[index].append(snp)
						gene_hits_counts[index].append(1)
						gene_hits_NonSyn_counts[index].append(1)
						gene_hits_Syn_counts[index].append(0)
						gene_hits_Amb_counts[index].append(0)
					if codon in gene_Codons[index]:
						codon_index = gene_Codons[index].index(codon)
						gene_Codons_counts[index][codon_index] += 1
						gene_Codons_NonSyn_counts[index][codon_index] += 1
					else:
						gene_Codons[index].append(codon)
						gene_Codons_counts[index].append(1)
						gene_Codons_NonSyn_counts[index].append(1)
						gene_Codons_Syn_counts[index].append(0)
						gene_Codons_Amb_counts[index].append(0)

				elif change_type == 'S':
					gene_counts_Syn[index] += 1
					if snp in gene_hits[index]:
						snp_index = gene_hits[index].index(snp)
						gene_hits_counts[index][snp_index] += 1
						gene_hits_Syn_counts[index][snp_index] += 1
					else:
						gene_hits[index].append(snp)
						gene_hits_counts[index].append(1)
						gene_hits_NonSyn_counts[index].append(0)
						gene_hits_Syn_counts[index].append(1)
						gene_hits_Amb_counts[index].append(0)
					if codon in gene_Codons[index]:
						codon_index = gene_Codons[index].index(codon)
						gene_Codons_counts[index][codon_index] += 1
						gene_Codons_Syn_counts[index][codon_index] += 1
					else:
						gene_Codons[index].append(codon)
						gene_Codons_counts[index].append(1)
						gene_Codons_NonSyn_counts[index].append(0)
						gene_Codons_Syn_counts[index].append(1)
						gene_Codons_Amb_counts[index].append(0)

				else:  #ambiguous
					gene_counts_Ambiguous[index] += 1
					if snp in gene_hits[index]:
						snp_index = gene_hits[index].index(snp)
						gene_hits_counts[index][snp_index] += 1
						gene_hits_Amb_counts[index][snp_index] += 1
					else:
						gene_hits[index].append(snp)
						gene_hits_counts[index].append(1)
						gene_hits_NonSyn_counts[index].append(0)
						gene_hits_Syn_counts[index].append(0)
						gene_hits_Amb_counts[index].append(1)
					if codon in gene_Codons[index]:
						codon_index = gene_Codons[index].index(codon)
						gene_Codons_counts[index][codon_index] += 1
						gene_Codons_Amb_counts[index][codon_index] += 1
					else:
						gene_Codons[index].append(snp)
						gene_Codons_counts[index].append(1)
						gene_Codons_NonSyn_counts[index].append(0)
						gene_Codons_Syn_counts[index].append(0)
						gene_Codons_Amb_counts[index].append(1)

			else:
				genes.append(gene)
				gene_counts.append(1)
				change_type = entry[14]
				if change_type == 'NS':
					gene_counts_NonSyn.append(1)
					gene_counts_Syn.append(0)
					gene_counts_Ambiguous.append(0)
					gene_hits.append([snp])
					gene_hits_counts.append([1])
					gene_hits_NonSyn_counts.append([1])
					gene_hits_Syn_counts.append([0])
					gene_hits_Amb_counts.append([0])
					gene_Codons.append([codon])
					gene_Codons_counts.append([1])
					gene_Codons_NonSyn_counts.append([1])
					gene_Codons_Syn_counts.append([0])
					gene_Codons_Amb_counts.append([0])

				elif change_type == 'S':
					gene_counts_NonSyn.append(0)
					gene_counts_Syn.append(1)
					gene_counts_Ambiguous.append(0)
					gene_hits.append([snp])
					gene_hits_counts.append([1])
					gene_hits_NonSyn_counts.append([0])
					gene_hits_Syn_counts.append([1])
					gene_hits_Amb_counts.append([0])
					gene_Codons.append([codon])
					gene_Codons_counts.append([1])
					gene_Codons_NonSyn_counts.append([0])
					gene_Codons_Syn_counts.append([1])
					gene_Codons_Amb_counts.append([0])

				else:  #ambiguous
					gene_counts_NonSyn.append(0)
					gene_counts_Syn.append(0)
					gene_counts_Ambiguous.append(1)
					gene_hits.append([snp])
					gene_hits_counts.append([1])
					gene_hits_NonSyn_counts.append([0])
					gene_hits_Syn_counts.append([0])
					gene_hits_Amb_counts.append([1])
					gene_Codons.append([codon])
					gene_Codons_counts.append([1])
					gene_Codons_NonSyn_counts.append([0])
					gene_Codons_Syn_counts.append([0])
					gene_Codons_Amb_counts.append([1])

	for i in range(len(genes)):
		output += genes[i] + '\t' + str(gene_counts[i]) + '\t' + str(gene_counts_Syn[i]) + '\t'
#		output += str(gene_counts_NonSyn[i]) + '\t' + str(gene_counts_Ambiguous[i]) + '\t' + str(len(gene_hits[i])) + '\t' + str(len(gene_Codons[i])) + '\t'
		output += str(gene_counts_NonSyn[i]) + '\t' + str(gene_counts_Ambiguous[i]) + '\t' + str(len(gene_hits[i])) + '\t' + str(len(gene_Codons[i])) + '\n'
#		gene_hits_out = ''
#		gene_hits_counts_out = ''
#		gene_hits_NonSyn_counts_out = ''
#		gene_hits_Syn_counts_out = ''
#		gene_hits_Amb_counts_out = ''
#		for j in range(len(gene_hits[i])):
#			gene_hits_out += str(gene_hits[i][j])
#			gene_hits_counts_out += str(gene_hits_counts[i][j])
#			gene_hits_NonSyn_counts_out += str(gene_hits_NonSyn_counts[i][j])
#			gene_hits_Syn_counts_out += str(gene_hits_Syn_counts[i][j])
#			gene_hits_Amb_counts_out += str(gene_hits_Amb_counts[i][j])
#			if j != len(gene_hits[i])-1:
#				gene_hits_out += ","
#				gene_hits_counts_out += ","
#				gene_hits_NonSyn_counts_out += ","
#				gene_hits_Syn_counts_out += ","
#				gene_hits_Amb_counts_out += ","
#		output += gene_hits_out +'\t'+ gene_hits_counts_out + '\t' + gene_hits_NonSyn_counts_out + '\t' + gene_hits_Syn_counts_out + '\t' + gene_hits_Amb_counts_out + '\t'
#		gene_Codons_out = ''
#		gene_Codons_counts_out = ''
#		gene_Codons_NonSyn_counts_out = ''
#		gene_Codons_Syn_counts_out = ''
#		gene_Codons_Amb_counts_out = ''
#		for j in range(len(gene_Codons[i])):
#			gene_Codons_out += str(gene_Codons[i][j])
#			gene_Codons_counts_out += str(gene_Codons_counts[i][j])
#			gene_Codons_NonSyn_counts_out += str(gene_Codons_NonSyn_counts[i][j])
#			gene_Codons_Syn_counts_out += str(gene_Codons_Syn_counts[i][j])
#			gene_Codons_Amb_counts_out += str(gene_Codons_Amb_counts[i][j])
#			if j != len(gene_Codons[i])-1:
#				gene_Codons_out += ","
#				gene_Codons_counts_out += ","
#				gene_Codons_NonSyn_counts_out += ","
#				gene_Codons_Syn_counts_out += ","
#				gene_Codons_Amb_counts_out += ","
#		output += gene_Codons_out + '\t' + gene_Codons_counts_out + '\t' + gene_Codons_NonSyn_counts_out + '\t' + gene_Codons_Syn_counts_out + '\t' + gene_Codons_Amb_counts_out + '\n'

	return output

def countMutationEvents(arguments):
#	output = "Gene_Tag\tME\tSyn\tPA\tAmb\tSites\tCodons\tSites_list\tSite_Counts\tPA_Site_Counts\tSyn_Site_Counts\tAmb_Site_Counts\tCodons_list\tCodon_Counts\tPA_Codon_Counts\tSyn_Codon_Counts\tAmb_Codon_Counts\n"
	output = "Gene_Tag\tME\tSyn\tPA\tAmb\tSites\tCodons\n"
	mutation_events_file = readInput(arguments.input_file)
	output = countEvents(mutation_events_file,output)
	writeOutput(arguments.output_prefix+"ME_counts_by_gene.txt",output)
	return

def main():
	countMutationEvents(parseArguments())
	return

if __name__ == '__main__':
	main()
