#!/bin/bash
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --output=./lbpWell.out
#SBATCH --error=./lbpWell.err

/home2/azaritsky/code/applications/2dActionRecognition/jobDistributionScriptLiveHistology/matlab/job_wrapper_lbp.sh /home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/ pcAccumulateLBPWell 1>/home2/azaritsky/logsBioHPC/LCH/log/lbpWell.std.out 2>/home2/azaritsky/logsBioHPC/LCH/log/lbpWell.std.err
