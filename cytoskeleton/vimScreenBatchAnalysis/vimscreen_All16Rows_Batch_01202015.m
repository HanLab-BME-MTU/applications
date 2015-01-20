#!/bin/bash
 
#SBATCH --job-name=Liya_test_RowWithInd
#SBATCH --partition=384GB
#SBATCH --nodes=1
#SBATCH --time=0
#SBATCH --mail-type=ALL
#SBATCH --output=myjob.%j.%N.out
#SBATCH --error=myjob.%j.%N.err 

# Load matlab 2013a using module
module add matlab/2013a

# Run matlab in batch mode
# run commands specified with -r option
# direct output info specified -logfile
# run as a background process (&)
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(1,[]), exit" -logfile "matlab_script_01202015_1_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(2,[]), exit" -logfile "matlab_script_01202015_2_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(3,[]), exit" -logfile "matlab_script_01202015_3_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(4,[]), exit" -logfile "matlab_script_01202015_4_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(5,[]), exit" -logfile "matlab_script_01202015_5_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(6,[]), exit" -logfile "matlab_script_01202015_6_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(7,[]), exit" -logfile "matlab_script_01202015_7_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(8,[]), exit" -logfile "matlab_script_01202015_8_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(9,[]), exit" -logfile "matlab_script_01202015_9_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(10,[]), exit" -logfile "matlab_script_01202015_10_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(11,[]), exit" -logfile "matlab_script_01202015_11_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(12,[]), exit" -logfile "matlab_script_01202015_12_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(13,[]), exit" -logfile "matlab_script_01202015_13_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(14,[]), exit" -logfile "matlab_script_01202015_14_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(15,[]), exit" -logfile "matlab_script_01202015_15_output.txt" &
matlab -singleCompThread -nodisplay -r "script_FilamentRun_Rows_withInd_SingleThread_test(16,[]), exit" -logfile "matlab_script_01202015_16_output.txt" &

# wait for all background processes to finish
wait

 
