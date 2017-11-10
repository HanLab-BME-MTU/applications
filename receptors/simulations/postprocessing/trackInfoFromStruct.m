function [xyCoordAll,xyCoordM,xyCoordS,trackStartRow] = trackInfoFromStruct(tracks)

%get number of tracks

numTracks = length(tracks);

%get number of frames
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%initialize output variables
xyCoordAll = [];
xyCoordM = [];
xyCoordS = [];
trackStartRow = zeros(numTracks,1);

%go over all tracks ...
for i = 1 : numTracks
    
    %get sequence of events of track
    seqOfEvents = tracks(i).seqOfEvents;
    
    %get start and end times of track
    startTime = seqOfEvents(1,1);
    endTime   = seqOfEvents(end,1);
    
    %extract track coordinates from structure
    tracksCoordAmpCG = tracks(i).tracksCoordAmpCG;
    xyCoordTmp = zeros(size(tracksCoordAmpCG,1),numFrames*2);
    xyCoordTmp(:,2*(startTime-1)+1:2:2*(endTime-1)+1) = tracksCoordAmpCG(:,1:8:end);
    xyCoordTmp(:,2*(startTime-1)+2:2:2*(endTime-1)+2) = tracksCoordAmpCG(:,2:8:end);
    
    %determine row where this compound track starts in xyCoordAll
    trackStartRow(i) = size(xyCoordAll,1) + 1;
    
    %get merge events
    mergeIndx = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));
    
    %go over all merges
    for iMerge = mergeIndx'
        
        %get merge time
        timeMerge = seqOfEvents(iMerge,1);
        
        %store coordinates belonging to this merge in xyCoordM
        xyCoordM = [xyCoordM; [timeMerge ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-1)+1:2*(timeMerge-1)+2) ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-2)+1:2*(timeMerge-2)+2) ...
            xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-2)+1:2*(timeMerge-2)+2)]]; %#ok<AGROW>
        
        %add coordinates after merge to merging track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-1)+1:2*(timeMerge-1)+2) = ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-1)+1:2*(timeMerge-1)+2);
        
    end
    
    %get split events
    splitIndx = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)));
    
    %go over all splits
    for iSplit = splitIndx'
        
        %get split time
        timeSplit = seqOfEvents(iSplit,1);
        
        %store coordinates belonging to this split in xyCoordS
        xyCoordS = [xyCoordS; [timeSplit ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2) ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-1)+1:2*(timeSplit-1)+2) ...
            xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-1)+1:2*(timeSplit-1)+2)]]; %#ok<AGROW>
        
        %add coordinates before split to splitting track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-2)+1:2*(timeSplit-2)+2) = ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2);
        
    end
    
    %store coordinates in xyCoordAll
    xyCoordAll = [xyCoordAll; xyCoordTmp]; %#ok<AGROW>
    
end

%replace NaNs by zero
xyCoordAll(isnan(xyCoordAll)) = 0;
