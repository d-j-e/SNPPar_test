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
python compareResults.py -e <expected> 
						 -o <observed>
						 -c <total snp count> 
						 -p <population> 
						 -s <sample size> 
						 -t <run time>
						 -m <memory use>
'''
#
# Last modified - 13/01/2020
#

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nCompares expected and observed homoplasic mutation events')
	parser.add_argument('-e', '--expected', type=str, required=True, help='Expected homoplasic mutation events - from simulation (required)')
	parser.add_argument('-o', '--observed', type=str, required=True, help='Observed homoplasic mutation events - from SNPPar (required)')
	parser.add_argument('-c', '--total_SNP_count', type=int, required=True, help='Total number of SNPs tested (required)')
	parser.add_argument('-p', '--population', type=str, required=True, help='Population tested eg. HCMC_L2 (required)')
	parser.add_argument('-s', '--size', type=int, required=True, help='Sample size (required)')
	parser.add_argument('-r', '--replicate', type=str, required=True, help='replicate number e.g. r1 (required)')
	parser.add_argument('-t', '--time', type=float, required=True, help='Time taken to run in seconds (required)')
	parser.add_argument('-m', '--memory', type=int, required=True, help='Memory use in bytes (required)')
	parser.add_argument('-T', '--tree', type=str, required=True, help='Phylogenetic tree (required)')
	parser.add_argument('-S', '--sorting', type=str, required=True, help='sorting used: simple (S), intermediate (I) or complex (C) (required)')
	return parser.parse_args()

def readInput(input_file_name):
	input_file_handle = open(input_file_name, 'r')
	input_file = input_file_handle.readlines()
	input_file_handle.close()
	return input_file

def outputResults(output_file_name,results):
	output_file_handle = open(output_file_name,"a")
	output = ''
	for i in range(len(results)):
		output += str(results[i])
		if i < len(results)-1:
			output += '\t'
	output +='\n'
	output_file_handle.write(output)
	output_file_handle.close()
	return

def readTree(new_tree):
	tree = Tree(new_tree,format=1)
	edge = 1
	for node in tree.traverse("preorder"):
		if not node.is_leaf():
			node.name = "N%d" %edge
			edge += 1
	return tree

def compareHomoplasies(test_set_expected,test_set_observed,total_SNP_count,population,size,time,memory,tree,output_file_name1,output_file_name2,sorting):
	test_tree = readTree(tree)
	#set time to minutes
	time = round(time/60,2)
	#set memory to megabytes
	memory = round(memory/1024/1024,2)
	# unique homoplasic snps from expected
	unique_expected = []
	for i in range(len(test_set_expected)):
		new = test_set_expected[i].split('\t')[0]
		if new not in unique_expected:
			unique_expected.append(new)
	# unique homoplasic snps from observed
	unique_observed = []
	for i in range(1,len(test_set_observed)):
		new = test_set_observed[i].split('\t')[0]
		if new not in unique_observed:
			unique_observed.append(new)
	True_Positive = 0
	False_Positive = 0
	False_Negative = 0
	unique_true_positive = []
	for item in unique_expected:
		if item in unique_observed: 
			True_Positive += 1
			unique_true_positive.append(item)
		else:
			False_Negative += 1
			outputResults(output_file_name1,['False Negative'])
			for i in range(len(test_set_expected)):
				if test_set_expected[i].split('\t')[0] == item:
					outputResults(output_file_name1,[test_set_expected[i].rstrip()])
	for item in unique_observed:
		if item not in unique_expected:
			False_Positive += 1
			outputResults(output_file_name1,['False Positive:'])
			for i in range(len(test_set_observed)):
				if test_set_observed[i].split('\t')[0] == item:
					outputResults(output_file_name1,[test_set_observed[i].rstrip()])
	True_Negative = total_SNP_count - (True_Positive + False_Positive + False_Negative)
	expected_type = []
	observed_type = []
	for item in unique_true_positive:
		SNP_test_set_expected = []
		for i in range(len(test_set_expected)):
			if test_set_expected[i].split('\t')[0] == item:
				SNP_test_set_expected.append(test_set_expected[i])
		SNP_test_set_observed = []
		for i in range(len(test_set_observed)):
			if test_set_observed[i].split('\t')[0] == item:
				SNP_test_set_observed.append(test_set_observed[i])
		expected_type = getHomoplasyType(SNP_test_set_expected,expected_type,True,test_tree)
		observed_type = getHomoplasyType(SNP_test_set_observed,observed_type,False,test_tree)
	type_result = compareHomoplasyType(expected_type,observed_type,output_file_name2)
	results = [sorting,population,size,total_SNP_count,time,memory,True_Positive,False_Positive,False_Negative,True_Negative] + type_result
	return results

def getHomoplasyType(SNP_test_set,ME_type,expected,tree):
	if expected:
		shift = 0
		for i in range(len(SNP_test_set)):
			SNP_test_set[i] = SNP_test_set[i].rstrip('\n')
	else:
		shift = 1
	types = ''
	root = 'N'
	for i in range(len(SNP_test_set)-1):
		for j in range(i+1,len(SNP_test_set)):
			#test for root node
			if SNP_test_set[i].split('\t')[1+shift] =='N1' or SNP_test_set[j].split('\t')[1+shift] =='N1':
				root = 'Y'
			#test pair for homoplasy
			if SNP_test_set[i].split('\t')[4+shift] == SNP_test_set[j].split('\t')[4+shift]:
				#parallel or convergent
				if SNP_test_set[i].split('\t')[3+shift] == SNP_test_set[j].split('\t')[3+shift]:
					#parallel
					if expected:
						types += 'P'
						#check not overlapping gene call
					elif SNP_test_set[i].split('\t')[6] == SNP_test_set[j].split('\t')[6]:
						types += 'P'
				else:
					#convergent
					types += 'C'
			elif (SNP_test_set[i].split('\t')[3+shift] == SNP_test_set[j].split('\t')[4+shift] 
				and SNP_test_set[i].split('\t')[4+shift] == SNP_test_set[j].split('\t')[3+shift]):
				if SNP_test_set[i].split('\t')[1+shift]== 'N1' or SNP_test_set[j].split('\t')[1+shift]== 'N1':
					if expected:
						types += 'R'
					elif SNP_test_set[i].split('\t')[6] == SNP_test_set[j].split('\t')[6]:
						types += 'R'
				else:
					node1 = tree.search_nodes(name=SNP_test_set[i].split('\t')[1+shift])[0]
					node2 = tree.search_nodes(name=SNP_test_set[j].split('\t')[1+shift])[0]
					if node1.get_common_ancestor(node2) == node1 or node2.get_common_ancestor(node1) == node2:
						if expected:
							types += 'R'
						elif SNP_test_set[i].split('\t')[6] == SNP_test_set[j].split('\t')[6]:
							types += 'R'
	ME_type.append([SNP_test_set[i].split('\t')[0],root,types])
	return ME_type

def compareHomoplasyType(expected_type,observed_type,output_file_name):
	parallel_events = 0
	convergent_events = 0
	revertant_events = 0
	correct_calls = 0
	incorrect_calls = 0
	root_node_correct = 0
	root_node_incorrect = 0
	single_calls = 0
	single_call_correct_P = 0
	single_call_correct_R = 0
	single_call_correct_C = 0
	single_call_incorrect_PR = 0
	single_call_incorrect_PC = 0
	single_call_incorrect_RP = 0
	single_call_incorrect_RC = 0
	single_call_incorrect_CP = 0
	single_call_incorrect_CR = 0
	for i in range(len(observed_type)):
		expected_type[i][2] = ''.join(sorted(expected_type[i][2]))
		observed_type[i][2] = ''.join(sorted(observed_type[i][2]))
		if expected_type[i][2] == observed_type[i][2]:
			correct_calls += 1
			parallel_events += expected_type[i][2].count('P')
			convergent_events += expected_type[i][2].count('C')
			revertant_events += expected_type[i][2].count('R')
			if expected_type[i][1] == 'Y' or observed_type[i][1] == 'Y':
				root_node_correct += 1
		else:
			if expected_type[i][1] == 'Y' or observed_type[i][1] == 'Y':
				root_node_incorrect += 1
			incorrect_calls += 1
			outputResults(output_file_name,['expected']+expected_type[i])
			outputResults(output_file_name,['observed']+observed_type[i])
		if len(expected_type[i][2]) == 1 and len(observed_type[i][2]) == 1:
			single_calls += 1
			if expected_type[i][2] == observed_type[i][2]:
				if expected_type[i][2] == 'P':
					single_call_correct_P += 1
				if expected_type[i][2] == 'R':
					single_call_correct_R += 1
				if expected_type[i][2] == 'C':
					single_call_correct_C += 1
			else:
				if expected_type[i][2] == 'P':
					if observed_type[i][2] == 'R':
						single_call_incorrect_PR += 1
					elif observed_type[i][2] == 'C':
						single_call_incorrect_PC += 1
				elif expected_type[i][2] == 'R':
					if observed_type[i][2] == 'P':
						single_call_incorrect_RP += 1
					elif observed_type[i][2] == 'C':
						single_call_incorrect_RC += 1
				if expected_type[i][2] == 'C':
					if observed_type[i][2] == 'P':
						single_call_incorrect_CP += 1
					elif observed_type[i][2] == 'R':
						single_call_incorrect_CR += 1
		else:
			print('SNP >2 ME: ' + expected_type[i][0])
			print(expected_type[i])
			print(observed_type[i])
	result = [correct_calls,incorrect_calls,parallel_events,convergent_events,revertant_events,root_node_correct,root_node_incorrect]
	result += [single_calls,single_call_correct_P,single_call_correct_C,single_call_correct_R,single_call_incorrect_PC,single_call_incorrect_PR]
	result += [single_call_incorrect_CP,single_call_incorrect_CR,single_call_incorrect_RP,single_call_incorrect_RC]
	return result

def testOutput(output_file_name,header):
	if not os.path.isfile(output_file_name):
		open(output_file_name, 'a').close()
		if header:
			outputResults(output_file_name,['Sorting','Population','Sample_Size','SNP_Calls','Time','Memory','True_Positive','False_Positive',
				'False_Negative','True_Negative','correct_calls','incorrect_calls','parallel_events','convergent_events','revertant_events',
				'root_node_correct','root_node_incorrect','single_calls','single_call_correct_P','single_call_correct_C','single_call_correct_R',
				'single_call_incorrect_PC','single_call_incorrect_PR','single_call_incorrect_CP','single_call_incorrect_CR',
				'single_call_incorrect_RP','single_call_incorrect_RC'])
	return

def compareResults(arguments):
	output_file_name = "homoplasic_test_results.tsv"
	output_file_name1 = arguments.sorting + '_' + arguments.population +'_'+ str(arguments.size) +'_'+ arguments.replicate + "_incorrect_homoplasic_calls.tsv"
	output_file_name2 = arguments.sorting + '_' + arguments.population +'_'+ str(arguments.size) +'_'+ arguments.replicate + "_incorrect_homoplasic_test_calls.tsv"
	testOutput(output_file_name,True)
	testOutput(output_file_name1,False)
	testOutput(output_file_name2,False)
	outputResults(output_file_name,compareHomoplasies(readInput(arguments.expected),readInput(arguments.observed),arguments.total_SNP_count,
		arguments.population,arguments.size,arguments.time,arguments.memory,arguments.tree,output_file_name1,output_file_name2,arguments.sorting))
	return

def main():
	compareResults(parseArguments())
	return

if __name__ == '__main__':
	main()
