#!/usr/bin/python2.7
#File name: pbs_script_dispatcher.py 
#NOTE: this file is supported in python 2.7
#NOTE: this file should be put in "LogNameExp"/job_submit,
#Discription: 
#	This script automatically submits some jobs to the queue and wait for a while, and check PBS queue to decide whether submit some jobs again or continue to wait.

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
total_nodes = 24 #number of physical nodes
wait_time = 300 #about half of a job finishing time is recommended.
start_job_index = 0
part = "64GB" # target partition
if len(sys.argv) == 1:
	print "starting job index = " + str(start_job_index) + " (by default)"
	print "waiting time = " + str(wait_time) + "s (by default)"
	print "target parition = " + part + " (by default)"
elif len(sys.argv) == 2:
	print "starting job index = " + sys.argv[1]
	print "waiting time = " + str(wait_time) + "s (by default)"
	print "target parition = " + part + " (by default)"
	wait_time = int(sys.argv[2])
elif len(sys.argv) == 3:
	print "starting job index = "  + sys.argv[1]
	print "waiting time = " + sys.argv[2]
	print "target parition = " + part + " (by default)"
	start_job_index = int(sys.argv[1])
	wait_time = int(sys.argv[2])
elif len(sys.argv) == 4:
	print "starting job index = "  + sys.argv[1]
	print "waiting time = " + sys.argv[2]
	print "target partition = " + sys.argv[3]
	start_job_index = int(sys.argv[1])
	wait_time = int(sys.argv[2])
	part = sys.argv[3]

#get a list of the script files
scripts = os.listdir('../autogen_script/')
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
	cmd = "sinfo | grep " + part +  " | grep idle | awk '{print $4}'"
	cmd_output = commands.getoutput(cmd)
	n_free_nodes = 0
	if cmd_output:
		n_free_nodes = int(cmd_output);
	print "#Current Free Nodes in " + part + " partition = " + str(n_free_nodes)
	if n_free_nodes == 0:
		print "Wait for more available nodes in " + part + " partition."
		time.sleep(wait_time)
		continue

	#submit jobs to queue until full	
	tmp_job_id = 0
	for i in range(0, n_free_nodes): #each bundle jobs=$total_nodes
		tmp_job_id = job_id + i
		if tmp_job_id >= total_num_scripts: #come to the end
			break
		print "SUBJOB_ID " +  str(tmp_job_id)
		script_name = scripts[tmp_job_id]
		cmd = 'sbatch ../autogen_script/' + script_name
		os.system(cmd) #dispatch a job
		print cmd
		cmd = 'echo ' + script_name +  ' >> submission.log'
		os.system(cmd) #keep a log for restarting from interruption.

	if tmp_job_id >= total_num_scripts: #come to the end
		break
	job_id = job_id + total_nodes;
	#after submission of a bunch of jobs, wait for a while to let PBS schedule the job to the nodes
	print "Wait for " + str(30) + " seconds to let SLURM react to the submission."
	time.sleep(30)

print str(num_scripts) + " job(s) submitted."
print "Submission finished."
#end here
