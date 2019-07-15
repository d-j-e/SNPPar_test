#!/usr/bin/env python
#
# 
# ****Function Description***
#
#
# Author(s) 
#			D. J. Edwards (David.Edwards@monash.edu) 
#
# Example command:
'''
python processSeqGen.py -i SeqGen_output -p prefix
'''
#
# Last modified - 
# Changes:
#

import os,sys,subprocess,string,re,random,collections,operator,argparse
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nIsolate List Randomiser ')
	parser.add_argument('-i', '--seqgen_output', type=str, required=True, help='SeqGen output file to process (required)')
	parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix for output files (required)')
	return parser.parse_args()

def isVariable(calls):
	test_call = calls[0]
	counter = 1
	while counter < len(calls):
		if not(test_call == calls[counter]):
			return True
		counter += 1
	return False

def tryIsInt(x):
	try:
		int(x)
	except:
		return False
	return True

def readSeqGenFile(seqgen_output):
	node_snp_calls = []
	node_names = []
	snps = []
	tip_snp_calls = []
	tip_names = []
	node_calls = []
	tip_calls = []
	input_file_handle = open(seqgen_output, 'r')
	line = input_file_handle.readline().lstrip()
	seq_count = int(line.split(' ')[0]) + 1
	seq_length = int(line.split(' ')[1].rstrip())
	tip_count = (seq_count + 1)//2
	for i in range(seq_length):
		node_calls.append([''])
		tip_calls.append([''])
	for i in range(seq_count):
		line_counter = 0
		line = input_file_handle.readline().rstrip()
		while not(line[line_counter]=='\t'):
			line_counter += 1
		name = line[0:line_counter]
		node_name = ''
		tip_name = ''
		if tryIsInt(name):
			if (int(name) >= tip_count + 1) and (int(name) <= seq_count):
				node_name = 'N' + str(int(name)-tip_count)
			else:
				tip_name = name
		else:
			tip_name = name
		if node_name:
			node_names.append(node_name)
			for j in range(seq_length):
				node_calls[j][0] = node_calls[j][0] + line[line_counter+1+j]
		else:
			tip_names.append(tip_name)
			for j in range(seq_length):
				tip_calls[j][0] = tip_calls[j][0] + line[line_counter+1+j]
	input_file_handle.close()
	for i in range(len(tip_calls)):
		if isVariable(tip_calls[i][0]):
			snps.append(i)
			tip_snp_calls.append(tip_calls[i])
			node_snp_calls.append(node_calls[i])
	return (node_snp_calls, node_names, snps, tip_snp_calls, tip_names)

def makeSNPtable(tip_names, snps, snp_calls):
	snptable = 'Pos'
	for name in tip_names:
		snptable += ',' + name 
	snptable += '\n'
	for i in range(len(snps)):
		snptable += str(snps[i]+1)
		for j in range(len(tip_names)):
			snptable += ',' + snp_calls[i][0][j] 
		snptable += '\n'
	return snptable

def tipOutput(prefix, snps, tip_snp_calls, tip_names):
	outputToFile(prefix+'.fasta',makeFasta(tip_snp_calls,tip_names))
	outputToFile(prefix+'_alleles.csv',makeSNPtable(tip_names, snps, tip_snp_calls))
	return

def makeFasta(snp_calls,names):
	fasta = ''
	for i in range(len(names)):
		fasta += '>' + names[i] + '\n'
		for j in range(len(snp_calls)):
			fasta += snp_calls[j][0][i]
		fasta += '\n'
	return fasta

def outputToFile(output_file_name, output):
	output_file_handle = open(output_file_name, "w")
	output_file_handle.write(output)
	output_file_handle.close()
	return

def nodeOutput(prefix, node_snp_calls, node_names):
	outputToFile(prefix+'_nodes.fasta',makeFasta(node_snp_calls,node_names))
	return

def output(prefix, node_snp_calls, node_names, snps, tip_snp_calls, tip_names):
	nodeOutput(prefix, node_snp_calls, node_names)
	tipOutput(prefix, snps, tip_snp_calls, tip_names)
	return

def processSeqGenFile(arguments):
	(node_snp_calls, node_names, snps, tip_snp_calls, tip_names) = readSeqGenFile(arguments.seqgen_output)
	output(arguments.prefix, node_snp_calls, node_names, snps, tip_snp_calls, tip_names)
	return

def main():
	processSeqGenFile(parseArguments())
	return

if __name__ == '__main__':
	main()
