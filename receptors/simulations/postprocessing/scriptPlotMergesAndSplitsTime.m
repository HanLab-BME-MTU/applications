%script to plot the mean times of merges and splits for simTracks and
%U-tracks

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160908/rD6/aP0p3/analysis/gap3/lR0p2/TrackingPackage/tracks';
saveDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160908/rD6/aP0p3/analysis/gap3/Results';
outDirNum = 2;

%load mean time values of splits for simulated data
                
        %iterate through the different runs
        for outDirIndx = 1 : length(outDirNum)
     tmp = load([sourceRoot,filesep,'msTimeInfoSim_10.mat']);
                timeMerge2SplitSim = tmp.msTimeInfoSim(outDirIndx).all.timeMerge2Split;
                clear tmp 
            %load mean time values of splits for U-track data
                 tmp = load([sourceRoot,filesep,'msTimeInfoUtra_10.mat']);
                timeMerge2SplitUtra = tmp.msTimeInfoUtra(outDirIndx).all.timeMerge2Split;
                clear tmp     
                    
                

subplot(2,1,1)       % add first plot in 2 x 1 grid
histogram(timeMerge2SplitSim);
title(['time Merge to Split Sim',int2str(outDirNum(outDirIndx))])

subplot(2,1,2)       % add second plot in 2 x 1 grid
histogram(timeMerge2SplitUtra);       % plot using + markers
title(['time Merge to Split U-track',int2str(outDirNum(outDirIndx))])

figH = gcf;
outFile = [saveDir,'/TimeMergeSplit_',int2str(outDirNum(outDirIndx))];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');

        end