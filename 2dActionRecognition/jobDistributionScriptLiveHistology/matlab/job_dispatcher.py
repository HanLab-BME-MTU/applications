#!/usr/bin/python2.7
#Name: job_dispatcher.py
#Author: Assaf Zaritsky (originally by Shenghua Wan)
#Date: 2015-04-23 (updated Dec. 2017)
#Description:
#	This script automatically submits jobs to the queue
#   NOTE: this file is supported in python 2.7


import os, re, sys, time
import commands #python 2

def tryint(s):
	try:
		return int(s)
	except:
		return s

def alphanum_key(s):
	""" Turn a string into a list of string and number chunks. 
		"z23a" -> ["z", 23, "a"]
	"""
	return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
	""" Sort the given list in the way that humans expect."""
	l.sort(key=alphanum_key)

#main function

#clean up previous log files, use CAUTION
print "Cleaning log folder ..."
filelist = [ f for f in os.listdir("/home2/azaritsky/logsBioHPC/LCH/log")  ]
for f in filelist:
    os.remove("/home2/azaritsky/logsBioHPC/LCH/log/" + f)

paramsList = []

# parse parameters file
print "Parsing parameters file ./params.txt ..."
with open("./params.txt", 'r') as f: # TODO: eventually move params.txt to a local location, so not updated in the repo..
	for line in f:
         line = line.rstrip()
         if not line.startswith('#') and len(line) > 0:
             paramsList.append(line)


partition = paramsList[4]
waitTime = int(paramsList[5])


start_job_index = 0

if len(sys.argv) == 1:
	print "starting job index = " + str(start_job_index) + " (by default)"
	print "waiting time = " + str(waitTime) + "s (by default)"
	print "target parition = " + partition + " (by default)"
elif len(sys.argv) == 2:
	print "starting job index = " + sys.argv[1]
	print "waiting time = " + str(waitTime) + "s (by default)"
	print "target parition = " + partition + " (by default)"
	waitTime = int(sys.argv[2])
elif len(sys.argv) == 3:
	print "starting job index = "  + sys.argv[1]
	print "waiting time = " + sys.argv[2]
	print "target parition = " + partition + " (by default)"
	start_job_index = int(sys.argv[1])
	waitTime = int(sys.argv[2])
elif len(sys.argv) == 4:
	print "starting job index = "  + sys.argv[1]
	print "waiting time = " + sys.argv[2]
	print "target partition = " + sys.argv[3]
	start_job_index = int(sys.argv[1])
	waitTime = int(sys.argv[2])
	partition = sys.argv[3]

#get a list of the script files
scripts = os.listdir('/home2/azaritsky/logsBioHPC/LCH/scripts/')
sort_nicely(scripts)
#print scripts
count = 0
total_num_scripts = len(scripts)
num_scripts = total_num_scripts - start_job_index 
print str(num_scripts) + " jobs to be submitted."
job_id = start_job_index

#NOTE
#submission policy: whenever the selected partition is available, use all the available nodes in it.
while job_id < total_num_scripts:
	#check free nodes
	cmd = "sinfo | grep " + partition +  " | grep idle | awk '{print $4}'"
	cmd_output = commands.getoutput(cmd)
	n_free_nodes = 0
	if cmd_output:
		n_free_nodes = int(cmd_output)
	print "#Current Free Nodes in " + partition + " partition = " + str(n_free_nodes)
	if n_free_nodes == 0:
		print "Wait for more available nodes in " + partition + " partition."
		time.sleep(waitTime)
		continue

	#submit jobs to queue until full	
	tmp_job_id = 0
	for i in range(0, n_free_nodes): 
		tmp_job_id = job_id + i
		if tmp_job_id >= total_num_scripts: #come to the end
			break
		print "SUBJOB_ID " +  str(tmp_job_id)
		script_name = scripts[tmp_job_id]
		cmd = 'sbatch /home2/azaritsky/logsBioHPC/LCH/scripts/' + script_name
		os.system(cmd) #dispatch a job
		print cmd
		cmd = 'echo ' + script_name +  ' >> /home2/azaritsky/logsBioHPC/LCH/log/submission.log'
		os.system(cmd) #keep a log for restarting from interruption.

	if tmp_job_id >= total_num_scripts: #come to the end
		break
	job_id = job_id + n_free_nodes
	#after submission of a bunch of jobs, wait for a while to let PBS schedule the job to the nodes
	print "Wait for " + str(30) + " seconds to let SLURM react to the submission."
	time.sleep(30)

print str(num_scripts) + " job(s) submitted."
print "Submission finished."
#end here
