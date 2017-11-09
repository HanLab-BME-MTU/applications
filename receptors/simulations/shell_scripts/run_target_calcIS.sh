#!/bin/sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --job-name=target_calcIS_sT25_dT0p01
#SBATCH --output=target_calcIS_sT25_dT0p01.out
#SBATCH --error=target_calcIS_sT25_dT0p01.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robel.yirdaw@utsouthwestern.edu


module load matlab


matlab -nodisplay -nosplash -noFigureWindows -logfile calcIS_rD4_aP0p5_lR0p2_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/targetISanalysis_sT25_dT0p01/; calcIS_cutOff('rD4/aP0p5/lR0p2',10,'target'); exit" &

matlab -nodisplay -nosplash -noFigureWindows -logfile calcIS_rD4_aP0p5_lR0p4_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/targetISanalysis_sT25_dT0p01/; calcIS_cutOff('rD4/aP0p5/lR0p4',10,'target'); exit" &

wait

