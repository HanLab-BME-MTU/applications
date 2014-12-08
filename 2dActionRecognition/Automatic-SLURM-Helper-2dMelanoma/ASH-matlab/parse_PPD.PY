#!/bin/python
#Name: 
#	parse_STPJD.py
#Date: 2014-02-04
#Description: 
#	this program reads file "Single_Thread_Program_Job_Description" and 
#generates SLURM job scripts for each node

import os, sys
import itertools #requires 2.6

#compatibility checking
if sys.version < '2.6':
	print 'Error: Python version > 2.6 is required.'
	print 'Current Python version = ' + sys.version
	sys.exit(1)

#clean up generated script, use CAUTION
filelist = [ f for f in os.listdir("./autogen_script")  ]
for f in filelist:
    os.remove("./autogen_script/" + f)

#process commandline arguemtns
part = '128GB'
if len(sys.argv) == 2:
	part = sys.argv[1]
else:
	print 'Use ' + part + ' partition (by default).'

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
				print "  line" + str(eff_line_no) + " :" + line #debug
				prog = line
				eff_line_no += 1
			elif(eff_line_no == 1):
				print "  line" + str(eff_line_no) + " :" + line #debug
				n_params = int(line)
				eff_line_no += 1	
			else:
				while eff_line_no >= 2 and eff_line_no < 2 + n_params:
					print "  line" + str(eff_line_no) + " :" + line #debug
					value_file_names.append(line)
					eff_line_no += 1
					break;
print
#read all value sets
set_list = []
for i in range(0, n_params):
	print "Reading finite value set: " + value_file_names[i]  + "..."
	value_set = []
	with open("./params/" + value_file_names[i]) as f:
		lines = f.readlines()
	value_set = [ x.rstrip() for x in lines ]
		#value_set = list(map(int, lines))
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
#Phase 2
#parse Single_Thread_Program_Job_Description
print "Parsing file ./Single_Thread_Program_Job_Description ..."
n_prog_per_node = 0
eff_line_no = 0
with open("./Single_Thread_Program_Job_Description", 'r') as f:
	for line in f:
		line = line.rstrip() #remove all trailing whitespaces
		if not line.startswith('#') and len(line) > 0:	
			if(eff_line_no == 0):
				print "  line " + str(eff_line_no) + " :" + line #debug
				n_prog_per_node = int(line)
				eff_line_no += 1
print '#progs per node = ', n_prog_per_node
#generate SLURM headers
slurm_header = '#!/bin/bash\n' \
			+ '#SBATCH --partition=128GB\n' \
			+ '#SBATCH --nodes=1\n' \
            + '#SBATCH --mail-user=Assaf.Zaritsky@utsouthwestern.edu\n' \
            + '#SBATCH --mail-type=all\n' #\
#			+ '#SBATCH --ntasks-per-node=' +  str(n_prog_per_node) + '\n'
			
			
if not os.path.exists('./autogen_script'):
	os.makedirs('./autogen_script')
#generate script for a node
n_whole_scripts = int(len(all_combs) / n_prog_per_node)
n_prog_last_script = len(all_combs) % n_prog_per_node
for i in range(0, n_whole_scripts):
	script_name = './autogen_script/job_' + str(i) + '.auto'
	print script_name
	with open(script_name, 'w') as f:
		sbatch_header = slurm_header + \
						'#SBATCH --output=auto_job_' + str(i) + '.out\n' + \
						'#SBATCH --error=auto_job_' + str(i) + '.err\n' \
						+ '\n'
		f.write(sbatch_header)
		for j in range(0, n_prog_per_node):
			comb_index = i * n_prog_per_node + j
			print comb_index
			#s =  prog + ' ' + ' '.join(map(str, all_combs[comb_index])) + ' & \n\n'
			log_name = 'task_' + '_'.join(map(str, all_combs[comb_index]))
			s =  prog + ' ' + ' '.join(map(str, all_combs[comb_index])) + ' 1>' + log_name + '.std.out' + ' 2>' + log_name + '.stderr'  + ' & \n\n'
			f.write(s)
		f.write('wait')
#last script
if(n_prog_last_script != 0):
	script_name = './autogen_script/job_' + str(n_whole_scripts) + '.auto'
	print script_name 
	with open(script_name, 'w') as f:
		sbatch_header = slurm_header + \
						'#SBATCH --output=auto_job_' + str(i) + '.out\n' + \
						'#SBATCH --error=auto_job_' + str(i) + '.err\n' \
						+ '\n'
		f.write(sbatch_header)
		for j in range(0, n_prog_last_script):
			comb_index = n_whole_scripts * n_prog_per_node + j
			print comb_index
			log_name = 'task_' + '_'.join(map(str, all_combs[comb_index]))
			s =  prog + ' ' + ' '.join(map(str, all_combs[comb_index])) + ' 1>' + log_name + '.stdout' + ' 2>' + log_name + '.stderr'  + ' & \n\n'
			f.write(s)
		f.write('wait')

print 'Summary'
print 'number of job scripts generated = ', n_whole_scripts, 'full script(s) and ',  \
		(1 if n_prog_last_script > 0 else 0), 'incomplete script(s)'
print 'each script contains up to ', n_prog_per_node, ' proceses'
