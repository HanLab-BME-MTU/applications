%Script to calculate interactions rates and cluster densities using
%aggregation information. Results will be saved in same directory as
%compound tracks and aggregation state.
%
%Khuloud Jaqaman, May 2015

sourceRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170425/varyDissRate/probe/probeIS_sT25_dT0p1';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170426/varyDissRate/probeISruns';

%Define strings for directory hierarchy as needed
  rDDir = {'rD14'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
 aPDir =  {'aP0p7'};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 dRDir={'dR1p75'};%,'dR1p5','dR1p75'
outDirNum =1:10;
lRDir = {'lR1'};%{'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'};
%define space and time information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: pay attention for the 'areaSideLen','timeStep',and 'sampleStep'
%values. For the low density simulation 'areaSideLen=25' and for high densities
%it is 12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% modification: Now infoSpaceTime has two new variables: systemState and
% 'intensityInfo'. systemState represents dynamic or static states.
% The values should be 1 for dynamic and 0 for static. Luciana de Oliveira, February 2017.

infoSpaceTime = struct('probDim',2,'areaSideLen',25,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10],'systemState',0);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    for dRDirIndx = 1 : length(dRDir)
    %iterate through the different runs
    for outDirIndx = 1 : length(outDirNum)
        for aPDirIndx = 1 : length(aPDir)
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,dRDir{dRDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),filesep,lRDir{lRDirIndx}];
                
                systemState=infoSpaceTime.systemState ;
                
                if systemState==1
                    %load compTracksAggregState
                    tmp = load([currDir,filesep,'compTracksAggregState.mat']);
                    compTracksAggregState = tmp.compTracksAggregState;
                elseif systemState==0
                    tmp = load([currDir,filesep,'detectionAggregState.mat']);
                    compTracksAggregState = tmp.detectionAggregState;
                end
                clear tmp
                
                %get rates and densities
                [rateOnPerClust,rateOffPerClust,densityPerClust,...
                    numClustForRateCalc,clustHistory,clustStats] = ...
                    clusterOnOffRatesAndDensityDynamicStatic(compTracksAggregState,infoSpaceTime);
                
                
                
                %save results
                
                saveDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,dRDir{dRDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
                
                save([saveDir,'/ratesAndDensity_dt0p1_T10'],'rateOnPerClust',...
                    'rateOffPerClust','densityPerClust','numClustForRateCalc',...
                    'clustHistory','clustStats','-v7.3');
                
                clear compTracksAggregState
                
            end %for each labelRatio
            
        end %for each outDir
        
    end %for each aP
    
    end
    
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx});
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
