function functionCalculateIntermidiateStatiscsMultipleParam(currFile,infoSpaceTime,cellIndex)

%This function calculates the intermediate statistics and the final output
%is in the format to calculate the Mahalanobis distance and p-value
%
% INPUT:
%
%
% currFile:         path where the cell array containing info about all the
%                   simulation and intermidiate statics parameters are saved.
%                   The cell array contains the folowing structure:
%                   column 1: file path;
%                   column 2: receptor density directory;
%                   column 3: dissociation rate directory;  
%                   column 4: association probability directory;                 
%                   column 5: labeled fraction directory;
%                   column 6: receptor density values;
%                   column 7: association probability values;
%                   column 8: dissociation rate values;
%                   column 9: labeled fraction values;
%                   column 10: output values;
% infoSpaceTime:    Structure with fields:
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
% Luciana de Oliveira, May 2017


%% load the cell array with all the information

% load the simulations info

temp=load([currFile,filesep,'cellInfoAllSim.mat']);
cellInfoAllSim=temp.cellInfoAllSim;
clear temp

%load the intermediate statistics info

temp=load([currFile,filesep,'cellInfoAllIntermidiateStatatistics.mat']);
cellInfoAllIntermidiateStatatistics=temp.cellInfoAllIntermidiateStatatistics;

destinationRoot =cellInfoAllIntermidiateStatatistics{1,1};

% for cellIndex
    
    %cellIndex=1:size(cellInfoAllSim,1)
    sourceRoot =cellInfoAllSim{cellIndex,1};
    
    %Define strings for directory hierarchy as needed
    
    rDDir =cellInfoAllSim{cellIndex,2};
    dRDir=cellInfoAllSim{cellIndex,3};
    aPDir =  cellInfoAllSim{cellIndex,4};
    lRDir = cellInfoAllSim{cellIndex,5};
    outDirNum =cellInfoAllSim{cellIndex,10};
    
    %% Calculations
    
    %step 1: calculate rates and density
    
    fprintf('\n Calculating rates and density. ')
    
    functionCalcMultipleRatesDensity(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)
    
%     %step 2: move rates and density
%     
     fprintf('\n Moving rates and density. ') 
%     
    functionMoveMultipleRatesDensity(sourceRoot,destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
%     
%     %step 3: collect rates and density
   fprintf('\n Collecting rates and density. ') 
%     
     functionCollectMultipleRatesDensity(destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
%     
%     %step 4: combine rates and density
%     
    fprintf('\n Combining rates and density. ') 
%     
      functionCombineMultipleRatesDensity(destinationRoot,rDDir,aPDir,dRDir,lRDir)
    
% end

