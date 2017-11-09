%Script to calculate interactions rates and cluster densities using
%aggregation information. Results will be saved in same directory as
%compound tracks and aggregation state.
%
%Khuloud Jaqaman, May 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/09/090514/probeISruns';
%Define strings for directory hierarchy as needed
rDDir = {'rD10'};
aPDir =  {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%,'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}
outDirNum =1:5;
lRDir = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5'};



%define space and time information
infoSpaceTime = struct('probDim',2,'areaSideLen',25,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10]);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
        
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
              
              saveDir=[sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                    mkdir(saveDir)
                    
                    
                save([saveDir,'/ratesAndDensity_dt0p1_T10'],'rateOnPerClust',...
                    'rateOffPerClust','densityPerClust','numClustForRateCalc',...
                    'clustHistory','clustStats','-v7.3');
                
                clear compTracksAggregState
                
            end %for each labelRatio
            
        end %for each outDir
        
    end %for each aP
    
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
