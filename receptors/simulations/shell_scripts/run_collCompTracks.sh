#!/bin/sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --job-name=collCT
#SBATCH --output=collCT.out
#SBATCH --error=collCT.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robel.yirdaw@utsouthwestern.edu

module load matlab

matlab -nodisplay -nosplash -noFigureWindows -logfile collCompTracks_log.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/targetISanalysis_sT25_dT0p01/; collectCompTracks; exit" &

wait
