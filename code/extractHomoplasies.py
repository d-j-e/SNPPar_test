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
python extractHomoplasies.py -h
'''
#
# Last modified - 08/01/2020
#

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nExtracts homoplasic mutation events (ME) from simulated data. Outputs list of ME at homoplasic SNPs')
	parser.add_argument('-n', '--nodes', type=str, required=True, help='nodes (FASTA) to process (required)')
	parser.add_argument('-t', '--tips', type=str, required=True, help='tips (FASTA) to process (required)')
	parser.add_argument('-T', '--tree_name', type=str, required=True, help='phylogenetic tree (Newick) (required)')
	parser.add_argument('-l', '--SNPList', type=str, required=True, help='List of SNP positions (required)')
	parser.add_argument('-p', '--output_prefix', type=str, default="", help='Prefix to add to output files')
	return parser.parse_args()

def readTree(tree_name):
	tree = Tree(tree_name)
	edge = 1
	node_names = []
	for node in tree.traverse("preorder"):
		if not node.is_leaf():
			node.name = "N%d" %edge
			node_names.append(node.name)
			edge += 1
	return (tree,node_names)

def readInput(input_file_name):
	input_file_handle = open(input_file_name, 'r')
	input_file = input_file_handle.readlines()
	input_file_handle.close()
	return input_file

def readFASTA(input_file_name):
	records = list(SeqIO.parse(input_file_name,"fasta"))
	return records

def outputMEList(output_file_name,ME_list):
	output_file_handle = open(output_file_name,"w")
	for item in ME_list:
		output = ''
		for i in range(len(item)):
			output += item[i]
			if i < len(item)-1:
				output += '\t'
		output += '\n' 
		output_file_handle.write(output)
	output_file_handle.close()
	return

def findSequence(name,records):
	for record in records:
		if record.id == name:
			return record.seq
	return

def getMutationEvents(tree,SNPList,tip_records,node_records,node_names):
	ME_list = []
	for node in tree.traverse("preorder"):
		if not node.is_leaf():
			node_name = node.name
			node_sequence = findSequence(node_name,node_records)
			for child in node.children:
				child_name = child.name
				if child_name in node_names: 
					child_sequence = findSequence(child_name,node_records)
					if not(child_sequence):
						print("node records")
						print(child_name)
				else:
					child_sequence = findSequence(child_name,tip_records)
					if not(child_sequence):
						print("tip records")
						print(child_name)
				for i in range(len(node_sequence)):
					if node_sequence[i] != child_sequence[i]:
						ME_list.append([SNPList[i].rstrip(),node_name,child_name,node_sequence[i],child_sequence[i]])
	for i in range(len(ME_list)):
		ME_list[i][0] = int(ME_list[i][0])
	ME_list = sorted(ME_list, key=itemgetter(0))
	for i in range(len(ME_list)):
		ME_list[i][0] = str(ME_list[i][0])
	return ME_list

def filterMutationEvents(ME_list,tree):
	last_item = -1
	homoplasy_list = []
	parallel_count = 0
	i = 0
	while i < len(ME_list):
		if not(last_item == -1) and ME_list[i][0] == ME_list[last_item][0]:		# same SNP position
			# find other mutation events with same position (sorted list!)
			same = 1
			while (i+same) < len(ME_list) and ME_list[i][0] == ME_list[i+same][0]:
				same += 1
			# test them pairwise
			homoplasic = False
			for j in range(last_item,last_item+same):
				for k in range(j+1,last_item+same+1):
					if ME_list[j][4] == ME_list[k][4]:		#same base change - either parallel or convergent
						homoplasic = True
					elif ME_list[j][3] == ME_list[k][4]:
						test_nodes = []
						node = tree.search_nodes(name=ME_list[j][1])[0]
						for item in node.get_descendants():
							test_nodes.append(item.name)
						for name in test_nodes:
							if name == ME_list[k][2]:
								homoplasic = True
					elif ME_list[j][4] == ME_list[k][3]:
						test_nodes = []
						node = tree.search_nodes(name=ME_list[k][1])[0]
						for item in node.get_descendants():
							test_nodes.append(item.name)
						for name in test_nodes:
							if name == ME_list[j][2]:
								homoplasic = True
			if homoplasic:
				for j in range(last_item,last_item+same+1):
					homoplasy_list.append(j)
			# reset i
			i = last_item + same
			last_item = i
			i += 1
		else: # reset
			last_item = i
			i += 1

	print("\nFound " + str(len(homoplasy_list)) + " mutation events involved with homoplasic SNPs\n")

	reduced_ME_list = []
	for i in homoplasy_list:
		reduced_ME_list.append(ME_list[i])
	return reduced_ME_list

def extract_hSNPs(arguments):
	homoplasic_output_file_name = arguments.output_prefix+"homoplasic_mutation_events.csv"
	(tree,node_names) = readTree(arguments.tree_name)
	SNPList = readInput(arguments.SNPList)
	tip_records = readFASTA(arguments.tips)
	node_records = readFASTA(arguments.nodes)
	outputMEList(homoplasic_output_file_name,filterMutationEvents(getMutationEvents(tree,SNPList,tip_records,node_records,node_names),tree))
	return

def main():
	extract_hSNPs(parseArguments())
	return

if __name__ == '__main__':
	main()
