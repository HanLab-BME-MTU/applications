#!/bin/bash
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --output=./lbp.out
#SBATCH --error=./lbp.err

/home2/azaritsky/code/applications/2dActionRecognition/jobDistributionScriptLiveHistology/matlab/job_wrapper_lbp.sh /home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/ pcAccumulateLBPNew pcMetaLBPNew 1>./lbp.std.out 2>./lbp.std.err
