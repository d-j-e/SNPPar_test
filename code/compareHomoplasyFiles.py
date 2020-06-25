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
python compareHomoplasyFiles.py -h
'''
#
# Last modified - 17/01/2020
#

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nCompares two lists of ME at homoplasic SNPs')
	parser.add_argument('-i', '--input_file_one', type=str, required=True, help='First ME file (required)')
	parser.add_argument('-I', '--input_file_two', type=str, required=True, help='Second ME file (required)')
	parser.add_argument('-n', '--name_file_one', type=str, required=True, help='Identifier for first file (required)')
	parser.add_argument('-N', '--name_file_two', type=str, required=True, help='Identifier for second file (required)')
	parser.add_argument('-p', '--output_prefix', type=str, default="", help='Prefix to add to output files')
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

def find_differences(arguments):
	output = ""
	file_one = readInput(arguments.input_file_one)
	file_two = readInput(arguments.input_file_two)
	name_file_one = arguments.name_file_one
	name_file_two = arguments.name_file_two
	i = 1
	j = 1
	while i < len(file_one) and j < len(file_two):
		if file_one[i] == file_two[j]:
			i += 1
			j += 1
		elif file_one[i].split('\t')[0] == file_one[j].split('\t')[0]:
			output += name_file_one + " and " +name_file_two+ " differ:\n"
			output += name_file_one + ":\n" + file_one[i]
			output += name_file_two + ":\n" + file_two[j]
			i += 1
			j += 1
		elif file_one[i].split('\t')[0] > file_one[j].split('\t')[0]:
			output += name_file_two + " differs:\n"
			output += file_two[j]
			j+=1
		else:
			output += name_file_one + " differs:\n"
			output += file_one[i]
			i+=1
	if i < len(file_one):
		output += name_file_one + " has extra call(s):\n"
		for k in range(i,len(file_one)):
			output += file_one[k]
	elif j < len(file_two):
		output += name_file_two + " has extra call(s):\n"
		for k in range(j,len(file_two)):
			output += file_two[k]
	writeOutput(arguments.output_prefix+"differences.txt",output)
	return

def main():
	find_differences(parseArguments())
	return

if __name__ == '__main__':
	main()
