function runReceptorAggregationSimple_HPC_PAR_(runIndex)
%RUNRECEPTORAGGREGATIONSIMPLE_HPC_PAR is used to perform simulations on the
%cluster. This is a template function and is called by a shell script. The
%function also handles calculation of rates.
%   
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
%   Robel Yirdaw, 10/31/14
% 

    fprintf('\n======================================');
    fprintf('\nRun # %d.',runIndex);
    fprintf('\n======================================\n');
    
    %Load random number seeds
    load('allRN_30.mat');

    modelParam = struct('diffCoef',0.1,'receptorDensity',10,'aggregationProb',[0.0;1.0;1.0;1.0;1.0;0.0;],...
        'aggregationDist',0.01,'dissociationRate',1.0,'labelRatio',1,'intensityQuantum',[1 0.3]);
    simParam = struct('probDim',2,'observeSideLen',25,'timeStep',0.01,'simTime',25,...
        'initTime',10,'randNumGenSeeds',allRN_30(runIndex));

    tic
    
    [receptorInfoAll,receptorInfoLabeled,~,~,assocStats,collProbStats] =...
        receptorAggregationSimple_new(modelParam,simParam);
         

    fprintf('\n====#1=====\n');
    currClock = clock(); 
    fprintf('\n%d-%d-%d %d:%d:%d\n',currClock(1),currClock(2),currClock(3),currClock(4),currClock(5),currClock(6));
    whos()
    fprintf('\n====#1=====\n');
    
    [compTracksALT,segmentStat] = ...
           aggregStateFromCompTracks_new(receptorInfoLabeled.compTracks,modelParam.intensityQuantum);  
       
    %clear compTracks
    
    %{
    [compTracksALT,segmentStat{simIndx}] = ...
            aggregStateFromCompTracks_noRepair(compTracks,modelParam.intensityQuantum); 
    %}    
    %Check aggregState against recept2clustAssign - the fourth return value
    %is compTracksR2C
    %{
    [mismatchInfo,aggregStateFromR2C,...
        aggregStateALL,~] =...
        checkAggregState(compTracksALT,recept2clustAssign,runIndex);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Begin rate calculations
    %NOTE: sending compTracksFromR2C not compTracks output from
    %aggregStateFromCompTracks.  FromR2C is mismatch free. (01/21/14)
    %clustHistoryAll = clusterHistoryFromCompTracks_aggregState(compTracksFromR2C);
    
    clustHistoryAll = clusterHistoryFromCompTracks_aggregState(compTracksALT.defaultFormatTracks);
    
    fprintf('\n====#2=====\n');
    %datestr(clock)
    currClock = clock(); 
    fprintf('\n%d-%d-%d %d:%d:%d\n',currClock(1),currClock(2),currClock(3),currClock(4),currClock(5),currClock(6));
    whos()
    fprintf('\n====#2=====\n');
    
    clear compTracksALT
    
    %Merge cluster histories into a single 2D array
    clustHistoryMerged = cat(1,clustHistoryAll{:,1});
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %Rates by cluster size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [rateAssocPerClust,rateDissocPerClust,eventTable,eventTable_mono] =...
        calcRatesByClustSize(clustHistoryMerged,simParam.timeStep);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Save association rate for monomers
    rateAssocMono = eventTable_mono(:,3);

    
    clear clustHistoryAll

    fprintf('\n====#3=====\n');
    %datestr(clock)
    currClock = clock(); 
    fprintf('\n%d-%d-%d %d:%d:%d\n',currClock(1),currClock(2),currClock(3),currClock(4),currClock(5),currClock(6));
    whos()
    fprintf('\n====#3=====\n');
    
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

