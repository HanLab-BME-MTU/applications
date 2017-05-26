function functionPreIntermidiateStatisticsMultipleParam(currFile,areaSideLen,intensityQuantum,cellIndex)

%function to calculate the pre-analysis for multiple parameters. In this
%function the file with all the information about the simulations will be
%load and calculate compTracks and aggregate states for all spectrum of
%parameters.
%
% INPUT
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
% areaSideLen:      Simulation/image side length values,
%                   which can be a single value or a value per
%                   side. In units of interest (e.g. um).
% intensityQuantum: Row vector with unit intensity mean
%                   and standard deviation (e.g. the intensity of a single
%                   fluorophore labeling a single receptor).
%
% Luciana de Oliveira, May 2017.

%% load the cell array with all the information

temp=load([currFile,filesep,'cellInfoAllSim.mat']);
cellInfoAllSim=temp.cellInfoAllSim;

%% Calculations

% for cellIndex=1:size(cellInfoAllSim,1)
    
    sourceRoot =cellInfoAllSim{cellIndex,1};
    
    %Define strings for directory hierarchy as needed
    
    rDDir =cellInfoAllSim{cellIndex,2};
    dRDir=cellInfoAllSim{cellIndex,3};
    aPDir =  cellInfoAllSim{cellIndex,4};
    outDirNum =cellInfoAllSim{cellIndex,10};
    lRDir = cellInfoAllSim{cellIndex,5};
    
    %step 0: check the simulations:
    %load the values of the parameters
    
    rDVal=cellInfoAllSim{cellIndex,6};
    aPVal=cellInfoAllSim{cellIndex,7};
    dRVal=cellInfoAllSim{cellIndex,8};
    
    % first test: if the input param are the same as output param
    
    fprintf('\n Checking the matching between inputs and outputs. ')
    
    simulationsWithProblems=functioncheckSimulations(sourceRoot,rDDir,aPDir,dRDir,outDirNum,rDVal,aPVal,dRVal);
    
    %if it pass in the first test, then go to the second test
    
    if cellfun(@isempty,simulationsWithProblems)
        
        fprintf('\n========================================================')
        
        %second test: if the density of receptors in the last frame is the same as in the input values.
        
        fprintf('\n Checking the matching between receptor densities.')
        
        simWithDensityProblems=functionCheckDensities(sourceRoot,rDDir,aPDir,dRDir,outDirNum,rDVal,areaSideLen);
    else
        
        error('There is a mismatch between input and output parameters.')
    end
    
    %if it pass in the second test, go to the compTracks calculation
    if cellfun(@isempty,simWithDensityProblems)
        
        %step 1: extract compTracks:
        
        fprintf('\n========================================================')
        
        fprintf('\n calculating compound tracks.')
        
%         functionExtractMultipleCompTracks(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
    else
        error('The density does not correspond to the input value')
    end
    %step 2: calculate aggregation state
    
    fprintf('\n========================================================')
    
    fprintf('\n calculating aggregate states from compTracks.')
    
%     functionCalcMultipleAggregState(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,intensityQuantum)
    
    
% end

