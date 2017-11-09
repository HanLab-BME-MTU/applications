#!/bin/sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --job-name=IS_sMin
#SBATCH --output=IS_sMin.out
#SBATCH --error=IS_sMin.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robel.yirdaw@utsouthwestern.edu

module load matlab


matlab -nodisplay -nosplash -noFigureWindows -logfile calcIS_sMin_rD6_aP0p3_lR0p2_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112514/targetIS_sT25_dT0p01_probeIS_sT25_dT0p01/; calcIS_sMin('rD4/aP0p5/lR0p2/'); exit" &

matlab -nodisplay -nosplash -noFigureWindows -logfile calcIS_sMin_rD6_aP0p3_lR0p4_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112514/targetIS_sT25_dT0p01_probeIS_sT25_dT0p01/; calcIS_sMin('rD4/aP0p5/lR0p4/'); exit" &


wait




