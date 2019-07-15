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
python compareResults.py -e <expected> -o <observed>  [-p <output_prefix>]
'''
#
# Last modified - 27/05/2019
# Changes: added parallel SNP filtering
#

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nExtracts mutation events (ME) from simulated data. Outputs include full list of ME and list of homoplasic ME')
	parser.add_argument('-e', '--expected', type=str, required=True, help='Expected parallel mutation events - from simulation (required)')
	parser.add_argument('-o', '--observed', type=str, required=True, help='Observed parallel mutation events - from "parallel" (required)')
	parser.add_argument('-t', '--total_SNP_count', type=int, required=True, help='Total number of SNPs tested (required)')
	parser.add_argument('-g', '--group', type=str, required=True, help='Group tested eg. HCMC_L2_r10p (required)')
	parser.add_argument('-r', '--replicate', type=int, required=True, help='Run replicate number (required)')
	parser.add_argument('-p', '--output_prefix', type=str, default="", help='Prefix to add to output files')
	parser.add_argument('-R', '--exclude_root', default=False, action="store_true" , help="Flag to exclude parallel SNPs with root as parent node (Default: False")
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

def getComparison(test_set_expected,test_set_observed,total_SNP_count,group,run):
	for item in test_set_expected:
		item[4] = sorted(item[4])
	for item in test_set_observed:
		item[4] = sorted(item[4])
	True_Positive = 0
	False_Positive = 0
	False_Negative = 0
	for item in test_set_expected:
		if item in test_set_observed: 
			True_Positive += 1
		else:
			False_Negative += 1
			if False_Negative == 1:
				print('False_Negative')
			print(item)
	for item in test_set_observed:
		if item not in test_set_expected:
			False_Positive += 1
			if False_Positive == 1:
				print('False_Positive')
			print(item)
	True_Negative = total_SNP_count - (True_Positive + False_Positive + False_Negative)
	Accuracy = (True_Positive+True_Negative)/total_SNP_count
	try:
		Specificity = True_Negative/(True_Negative+False_Positive)
	except:
		Specificity = None
	try:
		Sensitivity = True_Positive/(True_Positive+False_Negative)
	except:
		Sensitivity = None
	try:
		Precision = True_Positive/(True_Positive+False_Positive)
	except:
		Precision = None
	try:
		False_Positive_Rate = False_Positive/(True_Negative+False_Positive)
	except:
		False_Positive_Rate = None
	results = [run,group,True_Positive,False_Positive,False_Negative,True_Negative,Accuracy,Specificity,Sensitivity,Precision,False_Positive_Rate]
	return results

def in_test_set(setA,setB,test_set):
	set_A = setA.split('\t')
	set_B = setB.split('\t')	
	for item in test_set:
		if int(item[0]) == int(set_A[0]):
			if item[1] == set_A[3]:
				if item[2] == set_A[4].rstrip():
					if [set_A[1],set_A[2]] in item[4] and [set_B[1],set_B[2]] in item[4]:
						return True
	return False

def addToTestSet(setA,setB,test_set):
	addA = False
	addB = False
	index = None
	set_A = setA.split('\t')
	set_A[-1] = set_A[-1].rstrip() 
	set_B = setB.split('\t')
	set_B[-1] = set_B[-1].rstrip()	
	if test_set == []:
		test_set.append([set_A[0],set_A[3],set_A[4],1,[[set_A[1],set_A[2]],[set_B[1],set_B[2]]]])
	else:
		for i in range(len(test_set)):
			if int(test_set[i][0]) == int(set_A[0]):
				if test_set[i][1] == set_A[3]:
					if test_set[i][2] == set_A[4]:
						if [set_A[1],set_A[2]] not in test_set[i][4]:
							addA = True
						if [set_B[1],set_B[2]] not in test_set[i][4]:
							addB = True
							index = i
			else:
				addA = True
				addB = True
		if addA and addB:
			test_set.append([set_A[0],set_A[3],set_A[4],1,[[set_A[1],set_A[2]],[set_B[1],set_B[2]]]])
		elif addB:
			test_set[index][4].append([set_B[1],set_B[2]])
			test_set[index][3] += 1
	return test_set

def getTestSet(setA):
	test_set = []
	last_item = 0
	i = 1
	while i < len(setA):
		if setA[i].split('\t')[0] == setA[last_item].split('\t')[0]:
			same = 1
			while (i + same) < len(setA) and setA[i+same].split('\t')[0] == setA[last_item].split('\t')[0]:
				same += 1
			for j in range(last_item,last_item+same):
				for k in range(j+1,last_item+same+1):
					if setA[j].split('\t')[4] == setA[j].split('\t')[4]: #i.e. same base change
						if setA[j].split('\t')[3] == setA[k].split('\t')[3]: #i.e. same ancestor base
							if not(in_test_set(setA[j],setA[k],test_set)):
								test_set = addToTestSet(setA[j],setA[k],test_set)
			i = last_item + same
		last_item = i
		i += 1
	return test_set

def testOutput(output_file_name):
	if not os.path.isfile(output_file_name):
		open(output_file_name, 'a').close()
	return

def removeRootSNPs(SNP_count,expected,observed):
	expected_out = []
	observed_out = []
	snps_removed = []
	for item in expected:
		root_present = False
		for pair in item[4]:
			if pair[0] == 'N1': # i.e. root node
				root_present = True
		if not root_present:
			expected_out.append(item)
		else:
			snps_removed.append(item[0])
			SNP_count -= 1
	for item in observed:
		root_present = False
		for pair in item[4]:
			if pair[0] == 'N1': # i.e. root node
				root_present = True
		if not root_present:
			observed_out.append(item)
		else:
			if item[0] not in snps_removed:
				snps_removed.append(item[0])
				SNP_count -= 1
	return (SNP_count,expected_out,observed_out)

def compareResults(arguments):
	output_file_name = arguments.output_prefix+"parallel_test_results.tsv"
	testOutput(output_file_name)
	test_set_expected = getTestSet(readInput(arguments.expected))
	test_set_observed = getTestSet(readInput(arguments.observed))
	if arguments.exclude_root:
		(total_SNP_count, test_set_expected, test_set_observed) = removeRootSNPs(arguments.total_SNP_count, test_set_expected, test_set_observed)
	else:
		total_SNP_count = arguments.total_SNP_count

	outputResults(output_file_name,getComparison(test_set_expected,test_set_observed,total_SNP_count,arguments.group,arguments.replicate))
	return

def main():
	compareResults(parseArguments())
	return

if __name__ == '__main__':
	main()
