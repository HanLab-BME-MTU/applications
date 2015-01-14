#!/bin/bash
 
#SBATCH --job-name=Tony_test_RowWithInd
#SBATCH --partition=64GB
#SBATCH --nodes=1
#SBATCH --mem=16384 
#SBATCH --time=7-0:0:0
#SBATCH --mail-type=ALL
#SBATCH --output=myjob.%j.%N.out
#SBATCH --error=myjob.%j.%N.err
 
# Load matlab 2013a using module
module add matlab/2013a

# Run matlab in batch mode
# run commands specified with -r option
# direct output info specified -logfile
# run as a background process (&)
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(1,[]), exit" -logfile "matlab_script_01142015_1_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(2,[]), exit" -logfile "matlab_script_01142015_2_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(3,[]), exit" -logfile "matlab_script_01142015_3_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(4,[]), exit" -logfile "matlab_script_01142015_4_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(5,[]), exit" -logfile "matlab_script_01142015_5_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(6,[]), exit" -logfile "matlab_script_01142015_6_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(7,[]), exit" -logfile "matlab_script_01142015_7_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(8,[]), exit" -logfile "matlab_script_01142015_8_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(9,[]), exit" -logfile "matlab_script_01142015_9_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(10,[]), exit" -logfile "matlab_script_01142015_10_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(11,[]), exit" -logfile "matlab_script_01142015_11_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(12,[]), exit" -logfile "matlab_script_01142015_12_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(13,[]), exit" -logfile "matlab_script_01142015_13_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(14,[]), exit" -logfile "matlab_script_01142015_14_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(15,[]), exit" -logfile "matlab_script_01142015_15_output.txt" &
matlab -nodisplay -r "script_FilamentRun_Rows_withInd(16,[]), exit" -logfile "matlab_script_01142015_16_output.txt" &

# wait for all background processes to finish
wait

 
