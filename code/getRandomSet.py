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
python getRandomSet.py -l list -n number to return [-d output_directory]
'''
#
# Last modified - 
# Changes: 
#	15/07/2019	name change 
#

import os,sys,subprocess,string,re,random,collections,operator,argparse
from operator import itemgetter
from argparse import ArgumentParser

def parseArguments():
	parser = ArgumentParser(description='\nIsolate List Randomiser ')
	parser.add_argument('-l', '--isolate_list', type=str, required=True, help='Isolate list to sample (required)')
	parser.add_argument('-n', '--number', type=str, required=True, help='Number of isolates to return (required)')
	parser.add_argument('-o', '--ouput_name', type=str, required=True, help='Name of output file (required)')
	return parser.parse_args()

def readInput(input_file_name):
	input_file_handle = open(input_file_name, 'r')
	input_file = input_file_handle.readlines()
	input_file_handle.close()
	return input_file

def outputToFile(output_file_name, output):
	output_file_handle = open(output_file_name, "w")
	output_file_handle.write(output)
	output_file_handle.close()
	return

def main():
	# Parse arguments
	arguments = parseArguments()
	isolates = readInput(arguments.isolate_list)
	number = int(arguments.number)
	if number >= len(isolates):
		print("Cannot subsample " + str(number) + " from " + arguments.isolate_list)
	else:
		output = ""
		numbers = []
		for i in range(len(isolates)):
			numbers.append(i)
		random.shuffle(numbers)
		for i in range(number):
			output += isolates[numbers[i]]
		outputToFile(arguments.ouput_name, output)
	return

if __name__ == '__main__':
	main()
