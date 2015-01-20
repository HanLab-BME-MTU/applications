#!/bin/bash
 
#SBATCH --job-name=Liya_test_2RowWithInd
#SBATCH --partition=64GB
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=myjob.%j.%N.out
#SBATCH --error=myjob.%j.%N.err 

# Load matlab 2013a using module
module add matlab/2013a

# Run matlab in batch mode
# run commands specified with -r option
# direct output info specified -logfile
# run as a background process (&)
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(1,[]), exit" -logfile "matlab_script_0120B2015_1_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(2,[]), exit" -logfile "matlab_script_0120B2015_2_output.txt" &

# wait for all background processes to finish
wait

 
