%Script to calculate interactions rates and cluster densities using
%aggregation information. Results will be saved in same directory as
%compound tracks and aggregation state.
%
%Khuloud Jaqaman, May 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161028/target';

%Define strings for directory hierarchy as needed
rDDir = {'rD20'};
aPDir = {'aP0p5'};
outDirNum = 1:30;
lRDir = {'lR0p02','lR0p04'};


%define space and time information
infoSpaceTime = struct('probDim',2,'areaSideLen',25,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10]);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different runs
        for outDirIndx = 1 : length(outDirNum)
            
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
                %load compTracksAggregState
                tmp = load([currDir,filesep,'compTracksAggregState.mat']);
                compTracksAggregState = tmp.compTracksAggregState;
                clear tmp
                
                %get rates and densities
                [rateOnPerClust,rateOffPerClust,densityPerClust,...
                    numClustForRateCalc,clustHistory,clustStats] = ...
                    clusterOnOffRatesAndDensity(compTracksAggregState,infoSpaceTime);
                
                %save results
                save([currDir,'/ratesAndDensity_dt0p1_T10'],'rateOnPerClust',...
                    'rateOffPerClust','densityPerClust','numClustForRateCalc',...
                    'clustHistory','clustStats','-v7.3');
                
                clear compTracksAggregState
                
            end %for each labelRatio
            
        end %for each outDir
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
