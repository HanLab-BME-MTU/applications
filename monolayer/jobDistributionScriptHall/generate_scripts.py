#!/bin/python
#Name: generate_scripts.py
#Author: Assaf Zaritsky (originally by Shenghua Wan)
#Date: 2015-04-22
#Description: 
#	this program reads file "params.txt" and
#generates SLURM job scripts for each node

import os

#clean up generated script, use CAUTION
print "Cleaning pu scripts folder ..."
filelist = [ f for f in os.listdir("./scripts")  ]
for f in filelist:
    os.remove("./scripts/" + f)


paramsList = []

# parse parameters file
print "Parsing parameters file ./params.txt ..."
with open("./params.txt", 'r') as f:
	for line in f:
         line = line.rstrip()
         if not line.startswith('#') and len(line) > 0:
             paramsList.append(line)

path = paramsList[0]
userProg = paramsList[1]
nTasks = int(paramsList[2])
tasksPerNode = int(paramsList[3])
partition = paramsList[4]
waitTime = int(paramsList[5])
email = paramsList[6]

# parse user program
prog = userProg.split(".")[0]

if len(userProg.split(".")) > 1:
    ext = userProg.split(".")[1]
else:
    ext = '';

jobWrapperStr = path + '/job_wrapper.sh'
# progStr = path + '/' + userProg

if ext == 'm':
    s1 = jobWrapperStr + ' ' + path + ' ' + prog + ' '
else: # not debugged!
    s1 =  path + ' ' + prog + ' '

#generate SLURM headers
slurmHeader = '#!/bin/bash\n' \
			+ '#SBATCH --partition=' + partition + '\n' \
			+ '#SBATCH --nodes=1\n' \
            + '#SBATCH --mail-user=' + email + '\n' \
            + '#SBATCH --mail-type=all\n' #\
#			+ '#SBATCH --ntasks-per-node=' +  str(n_prog_per_node) + '\n'
			
			
if not os.path.exists('./scripts'):
	os.makedirs('./scripts')

#generate script for a node
nWholeScripts = int(nTasks / tasksPerNode)
nTasksLastScript = nTasks % tasksPerNode
for i in range(0, nWholeScripts):
	script_name = './scripts/job_' + str(i) + '.auto'
	print script_name
	with open(script_name, 'w') as f:
		sbatchHeader = slurmHeader + \
						'#SBATCH --output=./log/auto_job_' + str(i) + '_%j.out\n' + \
						'#SBATCH --error=./log/auto_job_' + str(i) + '_%j.err\n' \
						+ '\n'
		f.write(sbatchHeader)
		for j in range(0, tasksPerNode):
			comb_index = i * tasksPerNode + j
			print comb_index
			log_name = './log/task_' + str(i) + '_' + str(j+1) + '_' + str(comb_index+1)
			script =  s1 + str(comb_index+1) + ' 1>' + log_name + '.std.out' + ' 2>' + log_name + '.stderr'  + ' & \n\n'
			f.write(script)
		f.write('wait')

#last script
if nWholeScripts == 0:
    i = -1

if(nTasksLastScript != 0):
	script_name = './scripts/job_' + str(nWholeScripts) + '.auto'
	print script_name 
	with open(script_name, 'w') as f:
		sbatchHeader = slurmHeader + \
						'#SBATCH --output=./log/auto_job_' + str(nWholeScripts) + '_%j.out\n' + \
						'#SBATCH --error=./log/auto_job_' + str(nWholeScripts) + '_%j.err\n' \
						+ '\n'
		f.write(sbatchHeader)
		for j in range(0, nTasksLastScript):
			comb_index = nWholeScripts * tasksPerNode + j
			print comb_index
			log_name = './log/task_' + str(i+1) + '_' + str(j+1) + '_' + str(comb_index+1)
			script =  s1 + str(comb_index+1) + ' 1>' + log_name + '.std.out' + ' 2>' + log_name + '.stderr'  + ' & \n\n'
			f.write(script)
		f.write('wait')

print 'Summary'
print 'number of job scripts generated = ', nWholeScripts, 'full script(s) and ',  \
		(1 if nTasksLastScript > 0 else 0), 'incomplete script(s)'
print 'each script contains up to ', tasksPerNode, ' proceses'
