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
python SNPTableToList.py -i <SNPTable.csv> [-p <output_prefix>]
'''
#
# Last modified - 27/05/2019
# Changes:
#

import os,sys,subprocess,string,re,random,collections,operator,argparse
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nExtracts list of SNP positions from SNPTable ')
	parser.add_argument('-i', '--input', type=str, required=True, help='SNPTable to process (required)')
	parser.add_argument('-p', '--prefix', type=str, default="", help='Prefix to add to output files')
	return parser.parse_args()

def convertSNPTable(arguments):
	output_file_name = arguments.prefix + "SNPList.txt"
	f = open(arguments.input,'r')
	o = open(output_file_name, 'w')
	line = f.readline() #read header and discard...
	line = f.readline() #read first line of SNP calls
	while line:
		o.write(line.split(',')[0]+'\n')
		line = f.readline()
	f.close()
	o.close()
	return

def main():
	convertSNPTable(parseArguments())
	return

if __name__ == '__main__':
	main()
