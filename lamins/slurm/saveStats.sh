#!/bin/bash

# Created June 4th, 2015
# Mark Kittisopikul

#SBATCH --job-name saveStats
#SBATCH -p super
#SBATCH -N 1
#SBATCH -t 0-0:40:00
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user mark.kittisopikul@utsouthwestern.edu

module load matlab/2015a

matlab << MATLAB_CODE
disp('$1');
lamins.analysis.saveStats('$1');
exit
MATLAB_CODE
