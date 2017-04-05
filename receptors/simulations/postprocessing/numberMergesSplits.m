function [xyCoordAll0,numMerge0,numMerge1,numSplit0,numSplit1] = numberMergesSplits(tracksFinal,tracksSim)
%NUMBERMERGESSPLITS extracts xy coordinates of the simulations
% and the number of merges and splits in tracks %from Utrack and simulation
%
%SYNOPSIS  [numMerge0,numMerge1,numSplit0,numSplit1] = numberMergesSplits(tracksFinal,tracksSim)
% 
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : reformed tracks as obtained from the simulations.
%  
%       
%
%OUTPUT
%       xyCoordAll0 : x and y coordinates from the simulations.
%       numMerges0  : Number of merges from sim
%       numMerges1  : Number of merges from Utrack
%       numSplits0  : Number of splits from sim
%       numSplits0  : Number of splits from Utrack

%% process input variables

%extract track information out of tracksSim
[xyCoordAll0,xyCoordM0,xyCoordS0] = trackInfoFromStruct(tracksSim);

%extract track information out of tracksFinal
if isstruct(tracksFinal)
    
    [xyCoordAll1,xyCoordM1,xyCoordS1,~] = ...
        trackInfoFromStruct(tracksFinal);
    
       
else
    
    xyCoordAll1 = zeros(size(tracksFinal,1),size(tracksFinal,2)/4);
    xyCoordAll1(:,1:2:end) = tracksFinal(:,1:8:end);
    xyCoordAll1(:,2:2:end) = tracksFinal(:,2:8:end);
    xyCoordAll1(isnan(xyCoordAll1)) = 0;
   
    
    %do not look at merges and splits
    
    
end

%get number of frames in ground truth and in tracking results
numFrames0 = size(xyCoordAll0,2)/2;
numFrames1 = size(xyCoordAll1,2)/2;
if numFrames0 ~= numFrames1
 disp('ATTENTION: different number of frames in ground truth and simulation results!')
end

%convert coordinates in tracker output to coordinates in ground truth for
%later matching
for iFrame = 1 : numFrames1 
    
    % for merges
    if ~isempty(xyCoordM1)
        indxM1 = find(xyCoordM1(:,1)==iFrame);
        indxM2 = find(xyCoordM1(:,1)==iFrame+1);
        coord1 = [xyCoordM1(indxM1,2:3); xyCoordM1(indxM2,4:5); xyCoordM1(indxM2,6:7)];

            xyCoordM1(indxM1,2:3) = coord1(1:length(indxM1),:);
            xyCoordM1(indxM2,4:5) = coord1(length(indxM1)+1:length([indxM1;indxM2]),:);
            xyCoordM1(indxM2,6:7) = coord1(length([indxM1;indxM2])+1:end,:);
%         end
    end
    
    %repeat for splits
    if ~isempty(xyCoordS1)
        indxS1 = find(xyCoordS1(:,1)==iFrame);
        indxS2 = find(xyCoordS1(:,1)==iFrame+1);
        coord1 = [xyCoordS1(indxS1,4:5); xyCoordS1(indxS1,6:7); xyCoordS1(indxS2,2:3)];
           xyCoordS1(indxS1,4:5) = coord1(1:length(indxS1),:);
            xyCoordS1(indxS1,6:7) = coord1(length(indxS1)+1:length([indxS1;indxS1]),:);
            xyCoordS1(indxS2,2:3) = coord1(length([indxS1;indxS1])+1:end,:);
    end
    
end
%% check merges and splits

    
      
    %get number of merges in ground truth and tracking results
    numMerge0 = size(xyCoordM0,1);
    numMerge1 = size(xyCoordM1,1);

      
    %get number of splits in ground truth and tracking results
    numSplit0 = size(xyCoordS0,1);
    numSplit1 = size(xyCoordS1,1);
    
    
%% Subfunction 1

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

end
end
