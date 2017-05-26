
%function to calculate the pre-analysis for multiple parameters. In this
%function the file with all the information about the simulations will be
%load and calculate compTracks and aggregate states for all spectrum of
%parameters.
%
% INPUT
%
% currFile:       path where the cell array containing info about all the
%                 simulation parameters is saved. The cell array contains
%                 the folowing structure:
%                 column 1: file path;
%                 column 2: receptor density directory;
%                 column 3: association probability directory;
%                 column 4: dissociation rate directory;
%                 column 5: labeled fraction directory;
%                 column 6: receptor density values;
%                 column 7: association probability values;
%                 column 8: dissociation rate values;
%                 column 9: labeled fraction values;
%                 column 10: output values;

%load the cell array with all the information
temp=load([currFile,filesep,'cellInfoAllSimProbes.mat']);
cellInfoAllSimProbe=temp.cellInfoAllProbes;

%% Calculations

for cellIndex=1:size(cellInfoAllSimProbe,1)
    
    sourceRoot =cellInfoAllSimProbe{cellIndex,1};
    
    %Define strings for directory hierarchy as needed
    rDDir =cellInfoAllSimProbe{cellIndex,2};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
    aPDir =  cellInfoAllSimProbe{cellIndex,3};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
    dRDir=cellInfoAllSimProbe{cellIndex,4};%,'dR1p5','dR1p75'
    outDirNum =1:10;
    lRDir = cellInfoAllSimProbe{cellIndex,5};
    
    %step 0: check the simulations:
    %load the values of the parameters
    
    rDVal=cellInfoAllSimProbe{cellIndex,6};
    aPVal=cellInfoAllSimProbe{cellIndex,7};
    dRVal=cellInfoAllSimProbe{cellIndex,8};
    
    % first test: if the input param are the same as output param
    
    fprintf('\n Checking the matching between inputs and outputs.')
    
    simulationsWithProblems=functioncheckSimulations(sourceRoot,rDDir,aPDir,dRDir,outDirNum,rDVal,aPVal,dRVal);
    
    %if it pass in the first test, then go to the second test
    
    if cellfun(@isempty,simulationsWithProblems)
        
        fprintf('\n========================================================')
        
        %second test: if the density of receptors in the last frame is the same as in the input values.
        
        fprintf('\n Checking the matching between receptor density.')
        
        simWithDensityProblems=functionCheckDensities(sourceRoot,rDDir,aPDir,dRDir,outDirNum,rDVal,areaSideLen);
    else
        
        error('There is a mismatch between input and output parameters.')
    end
    
    %if it pass in the second test, go to the compTracks calculation
    if cellfun(@isempty,simWithDensityProblems)
        
        %step 1: extract compTracks:
        
        fprintf('\n========================================================')
        
        fprintf('\n calculating compound tracks.')
        
        functionExtractMultipleCompTracks(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
    else
        error('The density does not correspond to the input value')
    end
    %step 2: calculate aggregation state
    
    fprintf('\n========================================================')
    
    fprintf('\n calculating aggregate states from compTracks.')
    
    functionCalcMultipleAggregState(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,intensityQuantum)
    
    
end

