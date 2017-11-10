#!/bin/sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --job-name=probe_calcSingle_rD2
#SBATCH --output=probe_calcSingle_rD2.out
#SBATCH --error=probe_calcSingle_rD2.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robel.yirdaw@utsouthwestern.edu


module load matlab

aPDir=("aP0p2" "aP0p3" "aP0p4" "aP0p5" "aP0p6" "aP0p7" "aP0p8")
lRDir=("lR0p1" "lR0p2" "lR0p3" "lR0p4" "lR0p5" "lR0p6")

for i in {0..6}
do
   for j in {0..5}
   do
      for k in {1..5}
      do
         matlab -nodisplay -nosplash -noFigureWindows -logfile ${aPDir[$i]}/${lRDir[$j]}/calcSingle_log${k}.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/probeISanalysis_sT25_dT0p01/rD2/; calculateSingleRunQuants('${aPDir[$i]}/${lRDir[$j]}',$k); exit" &
      done
      wait

      for k in {6..10}
      do
         matlab -nodisplay -nosplash -noFigureWindows -logfile ${aPDir[$i]}/${lRDir[$j]}/calcSingle_log${k}.txt -r "cd /project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/probeISanalysis_sT25_dT0p01/rD2/; calculateSingleRunQuants('${aPDir[$i]}/${lRDir[$j]}',$k); exit" &
      done
      wait

   done
done

wait
