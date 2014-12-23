function calculateSingleRunQuants(fullPathToRuns,runIndex)
%CALCULATESINGLERUNQUANTS calculates and saves compTracksALT and rates for 
%each simulation. This function is part of the indirect inference based
%model calibration framework.
%
%  INPUT:
%           fullPathToRuns:     the location containing the repeated 
%                               simulation run outputs, saved as outx. 
%                               Each outx folder has compTracksx where x is
%                               the runIndex.
%           runIndex:           simulation number
%
%   OUTPUT:
%           compTracksALT are saved in a mat file as compTracksALTx.mat.
%           rates and other variables are saved in a file outx.mat
%
%
%  Robel Yirdaw, 09/11/14
%   

    fprintf('\n==============================================================================');
    fprintf('\nProcessing runs at %s and runIndex %d.',fullPathToRuns,runIndex);
    fprintf('\n==============================================================================');
    
    tic
        
    %Change to specified location
    cd([fullPathToRuns,'/out',int2str(runIndex)]);
    
    %Load compTracks
    loadStruct = load(['compTracks',int2str(runIndex),'.mat']);
    compTracks = loadStruct.compTracks;
    clear loadStruct
    
    %Track time and memory usage
    fprintf('\n====#1=====\n');
    currClock = clock(); 
    fprintf('\n%d-%d-%d %d:%d:%d\n',currClock(1),currClock(2),currClock(3),currClock(4),currClock(5),currClock(6));
    whos()
    fprintf('\n====#1=====\n');
    
    %Get ALT compTracks     
    [compTracksALT,segmentStat] = ...
           aggregStateFromCompTracks_new(compTracks,[1 0.3]);   
       
    clear compTracks
    
    %Can check aggregState against recept2clustAssign here. The fourth 
    %return value is compTracksR2C
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
    currClock = clock(); 
    fprintf('\n%d-%d-%d %d:%d:%d\n',currClock(1),currClock(2),currClock(3),currClock(4),currClock(5),currClock(6));
    whos()
    fprintf('\n====#2=====\n');
    
    
    %Merge cluster histories into a single 2D array
    clustHistoryMerged = cat(1,clustHistoryAll{:,1});
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %Rates by cluster size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [rateAssocPerClust,rateDissocPerClust,eventTable,eventTable_mono] =...
        calcRatesByClustSize(clustHistoryMerged,0.01);
    

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
    analysisETime = toc;
    
%%% 
    clear currClock
    
    %Save values to file.
    try       
        save(['compTracksALT',int2str(runIndex),'.mat'],'compTracksALT','-v7.3');
        clear compTracksALT
        
        outFile = ['out',int2str(runIndex),'.mat'];
        save(outFile,'-v7.3');
        fprintf('\n=====================================================');
        fprintf('\nVariables written to %s.',outFile);
        fprintf('\n=====================================================\n');
    catch outErr
        disp(outErr.message);
    end
        
    
end

