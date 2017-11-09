#!/bin/sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --job-name=IS_sp
#SBATCH --output=IS_sp.out
#SBATCH --error=IS_sp.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robel.yirdaw@utsouthwestern.edu

module load matlab


matlab -nodisplay -nosplash -noFigureWindows -logfile getIS_spVals_rD4_aP0p5_lR0p2_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112514/targetIS_sT25_dT0p01_probeIS_sT25_dT0p01/; getIS_spVals('rD4/aP0p5/lR0p2/'); exit" &

wait


matlab -nodisplay -nosplash -noFigureWindows -logfile getIS_spVals_rD4_aP0p5_lR0p2_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112514/targetIS_sT25_dT0p01_probeIS_sT25_dT0p01/; getIS_spVals('rD4/aP0p5/lR0p4/'); exit" &

wait



