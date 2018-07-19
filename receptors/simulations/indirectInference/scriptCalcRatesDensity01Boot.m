%Script to calculate interactions rates and cluster densities using
%aggregation information. Results will be saved in same directory as
%compound tracks and aggregation state.
%
%Khuloud Jaqaman, May 2015

sourceRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170327/varyDissRate/targetISruns';
bootDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170525/bootstrapping';
%Define strings for directory hierarchy as needed
rDDir = {'rD4'};%,'rD8','rD10','rD12','rD14','rD16'};
dRDir={'dR0p5'};
aPDir = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
outDirNum = 1:10;
lRDir ={'lR0p4'};%,'lR0p3','lR0p4','lR0p5'};
BootNumber=1:100;

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
                   dRDir{1},filesep, aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
                %load compTracksAggregState
                tmp = load([currDir,filesep,'compTracksAggregState.mat']);
                compTracksAggregState = tmp.compTracksAggregState;
                clear tmp
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% calculation clustHistory
       

 [clustHistoryAll,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState_new(compTracksAggregState.defaultFormatTracks);
                
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% bootstrapping
 
 %Number of bootstrapping repetitions
                for bootIndx=1:length(BootNumber)
                %get rates and densities
               [rateOnPerClust,rateOffPerClust,densityPerClust,...
    numClustForRateCalc,clustStats] = ...
    clusterOnOffRatesAndDensityBootstrapping(clustHistoryMerged,infoSpaceTime);

% create a directory to save the rates and densities
  currDirBoot = [bootDir,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                  mkdir(currDirBoot)

%save results
                save([currDirBoot,'/ratesAndDensity_dt0p1_T10_Boot_',int2str(BootNumber(bootIndx)),'.mat'],'rateOnPerClust',...
                    'rateOffPerClust','densityPerClust','numClustForRateCalc',...
                    'clustHistoryAll','clustStats','-v7.3');
                
               
                end %for each bootstrapping repetition
                
                 clear compTracksAggregState
                 
            end %for each labelRatio
            
        end %for each outDir
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
