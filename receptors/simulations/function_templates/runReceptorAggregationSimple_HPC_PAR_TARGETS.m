function runReceptorAggregationSimple_HPC_PAR_TARGETS(runIndex)
%RUNRECEPTORAGGREGATIONSIMPLE_HPC_PAR is used to perform target 
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
    fprintf('\nTarget IS run # %d.',runIndex);
    fprintf('\n======================================\n');
    
    %Load target random number seeds
    load('allRN_30_target.mat');

    modelParam = struct('diffCoef',0.1,'receptorDensity',4,'aggregationProb',[0.0;0.5;0.5;0.5;0.5;0.0;],...
        'aggregationDist',0.01,'dissociationRate',1.0,'labelRatio',[0.2;0.4;1.0],'intensityQuantum',[1 0.3]);
    simParam = struct('probDim',2,'observeSideLen',25,'timeStep',0.01,'simTime',25,...
        'initTime',10,'randNumGenSeeds',allRN_30_target(runIndex));

    tic
    
    [receptorInfoAll,receptorInfoLabeled,~,~,assocStats,collProbStats] =...
        receptorAggregationSimple_new(modelParam,simParam);
                 
    %070814
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

