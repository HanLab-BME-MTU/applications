#!/bin/python
#Name: 
#	parse_PPD.py
#Date: 2014-01-30
#Description: 
#	this program reads file "program_parameter_description" and 
#generates all combinations of possible parameter values.

import os, sys
import itertools #requires 2.6

#compatibility checking
if sys.version < '2.6':
	print 'Error: Python version > 2.6 is required.'
	print 'Current Python version = ' + sys.version
	sys.exit(1)
#main program
prog = ""	
n_params = ""
value_file_names = []
eff_line_no = 0 #effective/uncommented lines
#parse config file
print "Parsing file ./Program_Parameter_Description ..."
with open("./Program_Parameter_Description", 'r') as f:
	for line in f:
		line = line.rstrip() #remove all trailing whitespaces
		if not line.startswith('#') and len(line) > 0:
			
			if(eff_line_no == 0):
				print "  " + str(eff_line_no) + " :" + line #debug
				prog = line
				eff_line_no += 1
			elif(eff_line_no == 1):
				print "  " + str(eff_line_no) + " :" + line #debug
				n_params = int(line)
				eff_line_no += 1	
			else:
				while eff_line_no >= 2 and eff_line_no < 2 + n_params:
					print "  " + str(eff_line_no) + " :" + line #debug
					value_file_names.append(line)
					eff_line_no += 1
					break;
print
#read all value sets
set_list = []
n_param_values = [] #number of candidate parameter values
for i in range(0, n_params):
	print "Reading finite value set: " + value_file_names[i]  + "..."
	value_set = []
	with open("./params/" + value_file_names[i]) as f:
		lines = f.readlines()
	value_set = [ x.rstrip() for x in lines ]
	n_param_values.append(len(value_set))
	print value_set #debug
	set_list.append(value_set);
print "all value sets : "
print set_list #debug
#generate combinations: Cartesian product of all values sets
all_combs = set_list[0]
for i in range(1, n_params):
	comb = itertools.product(all_combs, set_list[i])
	all_combs = list(comb)
print all_combs #debug
#append the combinations to the binary
for comb in all_combs:
	print prog + ' ' + ' '.join(map(str, comb))

print 'Summary'
print 'number of parameters = ', n_params
print 'number of candidate parameter values', n_param_values
print 'number of parameter combinations = ', len(all_combs)