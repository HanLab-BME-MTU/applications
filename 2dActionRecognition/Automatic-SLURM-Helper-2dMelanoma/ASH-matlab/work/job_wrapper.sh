#!/bin/bash

module add matlab/2013a

cd /project/cellbiology/gdanuser/melanomaModel/Analysis/code/applications/2dActionRecognition/Automatic-SLURM-Helper-2dMelanoma/ASH-matlab/work

echo 'commandline arguments:'
#print arguments
for arg in $@
do
	echo $arg
done


#interpret the arguments according to the application interfaces
job_id=$1
task_id=$2

matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "runQLCHOneAtATimeMD($job_id,$task_id);"
