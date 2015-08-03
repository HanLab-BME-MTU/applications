#!/bin/bash

# Created June 4th, 2015
# Mark Kittisopikul

#SBATCH --job-name saveImages_20150503
#SBATCH -p 128GB
#SBATCH -N 1
#SBATCH -t 0-7:0:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user mark.kittisopikul@utsouthwestern.edu

module load matlab/2015a

matlab << MATLAB_CODE
disp('$1');
lamins.functions.saveImages('$1');
exit
MATLAB_CODE
