function [linkStats,fgapStats,bgapStats,gapType] = mtScoreLinksGaps(tracksFinal,tracksSim)
%SCORELINKSGAPSMS compares links, closed gaps and merges/splits obtained via tracking to the ground truth
%
%SYNOPSIS [linkStats,gapStats,mergeSplitStats] = scoreLinksGapsMS(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : Simulated tracks as obtained from
%                     simulateMimickCD36_MS (structure format).
%
%OUTPUT linkStats   : (Number of frames - 1) - by - 5 array. Row i 
%                     corresponds to the links from frame i to frame i+1.
%                     The columns show the number of:
%                     (1) ground truth links;
%                     (2) links resulting from tracking;
%                     (3) tracking links that are correct;
%                     (4) tracking links that are between real features but
%                         are wrong;
%                     (5) tracking links that involve a detection artifact
%                         (i.e. a detetion false positive).
%                     The sum of columns 3 through 5 = column 2.
%       fgapStats    : Two-column array with number of rows = number of
%                     forward gaps. First column shows gap length. Second
%                     column is
%                     0 - if gap is correctly closed;
%                     1 - if features just before and just after the gap
%                     are real but the connection is wrong;
%                     2 - if one or both features just before and just
%                     after the gap are detection artifacts.
%       bgapStats    : Two-column array with number of rows = number of
%                     forward gaps. First column shows gap length. Second
%                     column is
%                     0 - if gap is correctly closed;
%                     1 - if features just before and just after the gap
%                     are real but the connection is wrong;
%                     2 - if one or both features just before and just
%                     after the gap are detection artifacts.%       
% Sebastien Besson, June 2011
% Adapted from scoreLinksGaps.m

%% process input variables

%extract track information out of tracksSim
[xyCoordAll0,~,~,trackStartRow0] = trackInfoFromStruct(tracksSim);

if isstruct(tracksFinal)

    %extract track information out of tracksFinalnd merges/splits obtained via tracking to the ground truth
%
%SYNOPSIS [linkStats,gapStats,mergeSplitStats] = scoreLinksGapsMS(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matri
    [xyCoordAll1,~,~,trackStartRow1] = ...
        trackInfoFromStruct(tracksFinal);

else

    %directly extract track information out of tracksFinal
    xyCoordAll1 = zeros(size(tracksFinal,1),size(tracksFinal,2)/4);
    xyCoordAll1(:,1:2:end) = tracksFinal(:,1:8:end);
    xyCoordAll1(:,2:2:end) = tracksFinal(:,2:8:end);
    xyCoordAll1(isnan(xyCoordAll1)) = 0;
    trackStartRow1 = (1:size(xyCoordAll1,1))';

end

%get number of frames in ground truth and in tracking results
numFrames0 = size(xyCoordAll0,2)/2;
numFrames1 = size(xyCoordAll1,2)/2;
if numFrames0 ~= numFrames1
    numFrames0 = min(numFrames0,numFrames1);
    disp('ATTENTION: different number of frames in ground truth and simulation results!')
end

%% compare links (implicitly including merges and splits)

%initialize output variable
linkStats = NaN(numFrames0-1,5);
nd merges/splits obtained via tracking to the ground truth
%
%SYNOPSIS [linkStats,gapStats,mergeSplitStats] = scoreLinksGapsMS(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matri
%go over all frames ...
for iFrame = 1 : numFrames0 - 1
    
    %% number of links in ground truth - all correct by definition ...
    
    %get coordinates in ground truth
    tracksGT12 = xyCoordAll0(:,(iFrame-1)*2+1:(iFrame+1)*2);
    
    %get number of links in ground truth
    numLinksGT12 = length(find( tracksGT12(:,1) ~= 0 & tracksGT12(:,3) ~= 0 ));
    
    %% number of links from tracking code - some are correct, some are
    %% between real features but are wrong, others involve a false
    %% detection positive and thus are wrong ...
    
    %get coordinates in tracking results
    tracksCode12 = xyCoordAll1(:,(iFrame-1)*2+1:(iFrame+1)*2);
    
    %get number of links in tracking results
    numLinksTrackCode12 = length(find( tracksCode12(:,1)~=0 & tracksCode12(:,3)~=0 ));
            
    %% number of correct links from tracking code ...
    
    %get common links between ground truth and tracking results
    commonLinks = intersect(tracksGT12,tracksCode12,'rows');
    
    %remove zeros to get number of correct links
    numLinksCorrect12 = length(find( commonLinks(:,1)~=0 & commonLinks(:,3)~=0 ));
    
    %% number of wrong links from tracking code that are between real features ...
    
    %find the "real" features in the 2 frames
    realFeat1 = ismember(tracksCode12(:,1:2),tracksGT12(tracksGT12(:,1)~=0,1:2),'rows');
    realFeat2 = ismember(tracksCode12(:,3:4),tracksGT12(tracksGT12(:,3)~=0,3:4),'rows');
    
    %get number of all links from tracking code between real features
    numLinksReal12 = length(find( realFeat1~=0 & realFeat2~=0 ));
    
    %hence calculate number of wrong links between real features
    numLinksWrong12 = numLinksReal12 - numLinksCorrect12;
    
    %% number of links from tracking code involving false detection
    %% positives ...
        
    %get number of all links in tracking results that involve a false
    %positive (which are obviously wrong)
    numLinksFalse12 = length(find( ...
        (realFeat1 == 0 & tracksCode12(:,1) ~= 0 & tracksCode12(:,3) ~= 0) | ...
        (realFeat2 == 0 & tracksCode12(:,3) ~= 0 & tracksCode12(:,1) ~= 0) ));
    
    %% store numbers for output ...
    
    linkStats(iFrame,:) = [numLinksGT12 numLinksTrackCode12 numLinksCorrect12 ...
        numLinksWrong12 numLinksFalse12];
    
end

%% evaluate gaps in simulation

%find gaps in tracks from tracking code
gapInfo = findTrackGaps(tracksSim);
numGaps = size(gapInfo,1);

%allocate memory for gapStats
gapType = cell(numGaps,1);

%go over all gaps ...
for iGap = 1 : numGaps
    
    %get some gap information
    iTrack = gapInfo(iGap,1);
    iSegment = gapInfo(iGap,2);
    frameBef = gapInfo(iGap,3) - 1;
    gapLength = gapInfo(iGap,4);
    frameAft = frameBef + gapLength + 1;
    
    %find the coordinates from the tracking results before and after the gap
    xyCoordBef1 = xyCoordAll0(trackStartRow0(iTrack)+iSegment-1,(frameBef-1)*2+1:(frameBef-1)*2+2);
    xyCoordAft1 = xyCoordAll0(trackStartRow0(iTrack)+iSegment-1,(frameAft-1)*2+1:(frameAft-1)*2+2);
    displGap =xyCoordAft1-xyCoordBef1;
    
    %Determine if forward of backward gaps using the displacement of the
    %microtubule before the gap
    try
        xyCoordBef2 = xyCoordAll0(trackStartRow0(iTrack)+iSegment-1,(frameBef-2)*2+1:(frameBef-2)*2+2);
        displBef = xyCoordBef1-xyCoordBef2;
        gapOrientation = dot(displBef,displGap);
    catch ME1
        try
            xyCoordAft2 = xyCoordAll0(trackStartRow0(iTrack)+iSegment-1,(frameAft)*2+1:(frameAft)*2+2);
            displAft = xyCoordAft2-xyCoordAft1;
            gapOrientation = dot(displAft,displGap);
        catch ME2
            throw([ME1 ME2]);
        end
    end

    if gapOrientation>=0, 
        gapType{iGap} = 'f';
    else
        gapType{iGap} = 'b';
    end
end
   

%% evaluate gaps

%find gaps in tracks from tracking code
gapInfo = findTrackGaps(tracksFinal);
numGaps = size(gapInfo,1);

%allocate memory for gapStats
fgapStats = NaN(numGaps,2);
bgapStats = NaN(numGaps,2);

%go over all gaps ...
for iGap = 1 : numGaps
    
    %get some gap information
    iTrack = gapInfo(iGap,1);
    iSegment = gapInfo(iGap,2);
    frameBef = gapInfo(iGap,3) - 1;
    gapLength = gapInfo(iGap,4);
    frameAft = frameBef + gapLength + 1;
    
    %find the coordinates from the tracking results before and after the gap
    xyCoordBef1 = xyCoordAll1(trackStartRow1(iTrack)+iSegment-1,(frameBef-1)*2+1:(frameBef-1)*2+2);
    xyCoordAft1 = xyCoordAll1(trackStartRow1(iTrack)+iSegment-1,(frameAft-1)*2+1:(frameAft-1)*2+2);
    xyCoordBef2 = xyCoordAll1(trackStartRow1(iTrack)+iSegment-1,(frameBef-2)*2+1:(frameBef-2)*2+2);
    
    %determine whether either feature is a false positive
    realFeatBef = ismember(xyCoordBef1,xyCoordAll0(:,(frameBef-1)*2+1:(frameBef-1)*2+2),'rows');
    realFeatAft = ismember(xyCoordAft1,xyCoordAll0(:,(frameAft-1)*2+1:(frameAft-1)*2+2),'rows');
    
    if realFeatBef && realFeatAft
        
        %find featBef1 in frameBef in the ground truth tracks and look up which
        %feature it links to in the ground truth in frameAft
        indx = find( (xyCoordBef1(1) == xyCoordAll0(:,(frameBef-1)*2+1)) & ...
            (xyCoordBef1(2) == xyCoordAll0(:,(frameBef-1)*2+2)),1,'first' );
        
        xyCoordAft0 = xyCoordAll0(indx,(frameAft-1)*2+1:(frameAft-1)*2+2);

        %store wether gap is correctly closed (0) or not (1)
        gapStatsComponent = [gapLength ~min((xyCoordAft1 == xyCoordAft0))];

    else

        %indicate that gap involves detection artifacts
        gapStatsComponent = [gapLength 2];
        
    end
    
    %Determine if forward of backward gaps using the displacement of the
    %microtubule before the gap
    %Note: did not implement the case where length =1: should use the
    %displacement after gap as a safe mechanism
    displBef = xyCoordBef1-xyCoordBef2;
    displGap =xyCoordAft1-xyCoordBef1;
    if dot(displBef,displGap)>=0, 
        fgapStats(iGap,:) = gapStatsComponent;
    else
        bgapStats(iGap,:) = gapStatsComponent;
    end
end
   
% Remove Nan components in forward and backward gaps
fgapStats(isnan(fgapStats(:,1)),:)=[];
bgapStats(isnan(bgapStats(:,1)),:)=[];




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
            xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-2)+1:2*(timeMerge-2)+2)]];
        
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
            xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-1)+1:2*(timeSplit-1)+2)]];
        
        %add coordinates before split to splitting track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-2)+1:2*(timeSplit-2)+2) = ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2);

    end

    %store coordinates in xyCoordAll
    xyCoordAll = [xyCoordAll; xyCoordTmp];

end

%replace NaNs by zero
xyCoordAll(isnan(xyCoordAll)) = 0;
