#!/bin/bash

module add matlab/2013a

#cd /home2/azaritsky/HTC_scripts
cd /home2/azaritsky/HTC_scripts/slurm-helper-v2/work

echo 'commandline arguments:'
#print arguments
for arg in $@
do
	echo $arg
done


#interpret the arguments according to the application interfaces
job_id=$1
task_id=$2

matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "runQLCHOneAtATime($job_id,$task_id);"
