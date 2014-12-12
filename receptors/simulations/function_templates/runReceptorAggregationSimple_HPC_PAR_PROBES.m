function runReceptorAggregationSimple_HPC_PAR_PROBES(runIndex)
%RUNRECEPTORAGGREGATIONSIMPLE_HPC_PAR is used to perform probe 
%intermediate statistics simulations. This is a template function and is 
%called by a shell script.  
%
%   INPUT: 
%           runIndex:   the simulation (run) number
%
%   OUTPUT:
%           receptorInfoAll and receptorInfoLabeled are saved to separate
%           mat files.  All other quantities including rates are saved to
%           a single file, outx.mat, where x is runIndex. All are saved in
%           a directory call outx.
%
%  Robel Yirdaw, November 2014
% 
    fprintf('\n======================================');
    fprintf('\nProbe IS run # %d.',runIndex);
    fprintf('\n======================================\n');
    
    %Load probe random number seeds
    load('allRN_30.mat');

    modelParam = struct('diffCoef',0.1,'receptorDensity',2,'aggregationProb',[0.0;0.2;0.2;0.2;0.2;0.0;],...
        'aggregationDist',0.01,'dissociationRate',1.0,'labelRatio',[0.1;0.2;0.3;0.4;0.5;0.6;1.0],'intensityQuantum',[1 0.3]);
    simParam = struct('probDim',2,'observeSideLen',25,'timeStep',0.01,'simTime',25,...
        'initTime',10,'randNumGenSeeds',allRN_30(runIndex));

    fprintf('\n=========================================================================');
    fprintf('\ndiffCoeff = %g                           |   probDim = %g',modelParam.diffCoef,simParam.probDim);
    fprintf('\nreceptorDensity = %g                       |   observeSideLen = %g',modelParam.receptorDensity,simParam.observeSideLen);
    fprintf('\naggregationProb = [%g;%g;%g;%g;%g;%g]   |   timeStep = %g',modelParam.aggregationProb,simParam.timeStep);
    fprintf('\naggregationDist = %g                    |   simTime = %g',modelParam.aggregationDist,simParam.simTime);
    fprintf('\ndissociationRate = %g                      |   initTime = %g         ',modelParam.dissociationRate,simParam.initTime);
    fprintf('\nlabelRatio = [%g;%g;%g;%g;%g;%g;%g]  |   randNumGenSeed = %d   ',modelParam.labelRatio,simParam.randNumGenSeeds);
    fprintf('\nintensityQuantum = [%g %g]                |                         ',modelParam.intensityQuantum);
    fprintf('\n=========================================================================\n');
    
    tic
    
    [receptorInfoAll,receptorInfoLabeled,~,~,assocStats,collProbStats] =...
        receptorAggregationSimple_new(modelParam,simParam);
         
    elapsedTime = toc;
    
%%%
    try
        outDir = ['out',int2str(runIndex)];
        outFile = [outDir,'/out',int2str(runIndex),'.mat'];
        mkdir(outDir);
        
        %Save receptorInfoAll and receptorInfoLabeled separately
        save([outDir,'/receptorInfoAll',int2str(runIndex),'.mat'],'receptorInfoAll','-v7.3');
        clear receptorInfoAll

        save([outDir,'/receptorInfoLabeled',int2str(runIndex),'.mat'],'receptorInfoLabeled','-v7.3');
        clear receptorInfoLabeled
        
        save(outFile,'-v7.3');
        fprintf('\n=====================================================');
        fprintf('\nVariables written to %s.',outFile);
        fprintf('\n=====================================================\n');
    catch outErr
        disp(outErr.message);
    end
        
    
end

