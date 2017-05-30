
function functionCalcMultipleRatesDensityReceptorInfoAllStatic(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime,bootRepSize)
%function to calculate interactions rates and cluster densities using
%aggregation information. Results will be saved in same directory as
%compound tracks and aggregation state.
%SYNOPSIS  functionCalcMultipleRatesDensityReceptorInfoAllStatic(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)
%
%INPUT
%  sourceRoot : path for the receptorInfoLabeled.
%  rDDir: multiple receptor density directories that will be analyzed. It
%  is a cell in the form: {'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'}
%
%  aPDir      : association probability directories. A cell in the form: {'aP0p2','aP0p3','aP0p4'}
%
%  dRDir      : dissociation rate directories. A cell in the form: {'aP0p2','aP0p3','aP0p4'}
%
%  outDirNum  : list with the simulation number. 1:10 for the low density and
%  1:30 for the high density.
%
%  lRDir      : labed receptor directory. A cell in the form: {'lR0p2';'lR0p3';'lR0p4'}
%infoSpaceTime: Structure with fields:
%           .probDim        : Problem dimensionality.
%           .areaSideLen    : Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%           .timeStep       : Time between frames/time points. In units of
%                             interest (e.g. s).
%           .sampleStep     : Sampling time step, in same units as
%                             timeStep. Mostly relevant for
%                             simulated data where simulation time step
%                             might be 0.01 s but sampling time step of
%                             interest is e.g. 0.1 s.
%                             Optional. If not input, then sampleStep =
%                             timeStep.
%           .firstLastTP    : Row vector of first and last time points to
%                             use for calculating rates and densities. In
%                             same units as timeStep.
%                             If only one value is input, it is taken as
%                             the last time point.
%                             If no value is input, then all time points
%                             are used.
%
%
%            .systemState    : value saying if the analysis is for
%                              dynaminic or static data. Should be one of
%                              the two values:
%                              1-dynamic;
%                              0-static;
%OUTPUT
% There is no output, the results will be saved in the sourceRoot directory
% as /ratesAndDensity_dt0p1_T10.m
%
% Luciana de Oliveira, April 2017


fprintf('\n===============================================================');

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
                        
                        
%%                       %%%%%%%%%%%%%%%%%%%bootstrapping%%%%%%%%%%%%%%%%%%%%% 
                        for iboot=1: bootRepSize
                            compTracksAggregState = tmp.detectionAggregState{iboot};
                        
                                              
                        %get rates and densities
                        [rateOnPerClust{iboot},rateOffPerClust{iboot},densityPerClust{iboot},...
                            numClustForRateCalc{iboot},clustHistory{iboot},clustStats{iboot}] = ...
                            clusterOnOffRatesAndDensityDynamicStatic(compTracksAggregState,infoSpaceTime);
                        
                        end
                    end
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