#!/bin/sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --job-name=target_calcSingle_rD4
#SBATCH --output=target_calcSingle_rD4.out
#SBATCH --error=target_calcSingle_rDr.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robel.yirdaw@utsouthwestern.edu


module load matlab

aPDir=("aP0p5")
lRDir=("lR0p2" "lR0p4")

for i in {0..0}
do
   for j in {0..1}
   do
      for k in {1..5}
      do
         matlab -nodisplay -nosplash -noFigureWindows -logfile rD4/${aPDir[$i]}/${lRDir[$j]}/calcSingle_log${k}.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/targetISanalysis_sT25_dT0p01/; calculateSingleRunQuants('rD4/${aPDir[$i]}/${lRDir[$j]}',$k); exit" &
      done
      wait

      for k in {6..10}
      do
         matlab -nodisplay -nosplash -noFigureWindows -logfile rD4/${aPDir[$i]}/${lRDir[$j]}/calcSingle_log${k}.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/targetISanalysis_sT25_dT0p01/; calculateSingleRunQuants('rD4/${aPDir[$i]}/${lRDir[$j]}',$k); exit" &
      done
      wait

   done
done

wait
