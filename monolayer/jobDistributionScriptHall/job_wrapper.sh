#!/bin/bash

module add matlab/2015a

#interpret the arguments according to the application interfaces
directory=$1
prog=$2
ind=$3

# change to working directory (e.g, /home2/azaritsky/code/jobDistributionScript)
echo cd $directory
echo
cd $directory
echo

echo 'commandline arguments:'
#print arguments
for arg in $@
do
	echo $arg
done
echo

echo matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "$prog($ind), exit"
echo
matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "$prog($ind), exit"
