#!/usr/bin/python2.7
#Name: job_dispatcher.py
#Author: Assaf Zaritsky (originally by Shenghua Wan)
#Date: 2015-04-23
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
filelist = [ f for f in os.listdir("./log")  ]
for f in filelist:
    os.remove("./log/" + f)

paramsList = []

# parse parameters file
print "Parsing parameters file ./params.txt ..."
with open("./params.txt", 'r') as f:
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
scripts = os.listdir('./scripts/')
sort_nicely(scripts)
#print scripts
count = 0
total_num_scripts = len(scripts)
num_scripts = total_num_scripts - start_job_index 
print str(num_scripts) + " jobs to be submitted."
job_id = start_job_index

#NOTE
#submission policy:  submit all jobs, allow scheduler to manage jobs-submission
while job_id < total_num_scripts:
	#submit jobs to queue until full	
	print "SUBJOB_ID " +  str(job_id)
	script_name = scripts[job_id]
	cmd = 'sbatch ./scripts/' + script_name
	os.system(cmd) #dispatch a job
	print cmd
	cmd = 'echo ' + script_name +  ' >> ./log/submission.log'
	os.system(cmd) #keep a log for restarting from interruption.
	job_id = job_id + 1

print str(num_scripts) + " job(s) submitted."
print "Submission finished."
#end here
