#!/bin/bash
#SBATCH --job-name=serialJob                              # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-03:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=andrew.jamieson@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

module add matlab/2017a                                  # load software package

matlab -nodisplay -nodesktop -nosplash -r "run('/home2/s170480/matlab/applications/2dActionRecognition/explore/scripts/script_CellMD_OMETIFF_Gen2n3May15_arj.m'), exit"