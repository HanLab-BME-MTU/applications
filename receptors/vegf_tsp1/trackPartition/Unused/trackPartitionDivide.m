function tracksDivided = trackPartitionDivide(tracksPart,minSegLength)
%TRACKPARTITIONDIVIDE Takes the output structure of trackPartition and uses
%the partitioning information to cut up each track into segments each time
%the particle enters or exits a partition mask 
%   
%   tracksDivided = trackPartitionDivide(tracksPart,minSegLength)
%
%   Inputs:
%       tracksPart:     Output from trackPartition
%       minSegLength:   If more than 20% of segments are shorter than this
%                       length, this function will display a message.
%
%   Outputs:
%       tracksDivided:  Output structure with the following fields:
%       tracksFeatIndxCG, tracksCoordAmpCG:
%                       Same as in the trackCloseGapsKalmanSparse output, 
%                       however each partition segment is now on its own 
%                       line (as a new sub-track)
%       seqOfEvents:    Same as in the trackCloseGapsKalmanSparse output, 
%                       except now with partition "enter" and "exit" events 
%                       recorded. The first column indicates the frame of 
%                       the event, the second column can be
%                       1: Start of track segment
%                       2: End of track segment
%                       3: Particle entering partition mask (start of
%                       "inside" segment)
%                       4: Particle exiting partition mask (start of
%                       "outside" segment)
%                       The third column indicates which track segment is
%                       involved in the event, and the fourth column can
%                       be:
%                       NaN = start is a birth and end is a death,
%                       number = start is due to a split, end is due to a 
%                       merge, number is the index of track segment for the 
%                       merge/split.



nCompoundTracks = size(tracksPart,1);

% Output struct
tracksDivided = struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],...
    'seqOfEvents',[],'originSegment',[],'isInside',[]);
iRow = 1;

smallSegs = 0; % Keep track of the number of small segments
totalSegs = 0;
for i = 1:nCompoundTracks
    tracksIndxOld = tracksPart(i).tracksFeatIndxCG;
    tracksCoordOld = tracksPart(i).tracksCoordAmpCG;
    seqOfEventsOld = tracksPart(i).seqOfEvents;
    % Need to figure out how to preserve old track segment numbering
    nTracks = size(tracksIndxOld,1);
    seqPartition = tracksPart(i).seqPartition;
   
    tracksIndxNew = [];
    tracksCoordNew = [];
    seqOfEventsNew = seqOfEventsOld;
    isInsideNew = [];
    iTrack = 1;
        
    
    for j = 1:nTracks
            
        inside = tracksPart(i).inside{j};
        segmentsLeft = true;
        segStart = 1; % Start with at beginning of track

        
        while segmentsLeft
            if inside(segStart) == false % Current segment is outside
                segEnd = find(inside((segStart+1):end) == true,1,'first')+segStart;
                if isempty(segEnd) % If the end of the track is reached
                    segEnd = numel(inside);
                    segmentsLeft = false;
                end

                seqOfEventsTemp = [segStart,4,iTrack,NaN];
                seqOfEventsNew = [seqOfEventsNew; seqOfEventsTemp];
                isInsideNew = [isInsideNew; false];
            else % Current segment is inside
                segEnd = find(inside((segStart+1):end) == false,1,'first')+segStart;
                if isempty(segEnd) % If the end of the track is reached
                    segEnd = numel(inside);
                    segmentsLeft = false;
                end

                seqOfEventsTemp = [segStart,3,iTrack,NaN];
                seqOfEventsNew = [seqOfEventsNew; seqOfEventsTemp];
                isInsideNew = [isInsideNew; true];
            end
            if (segEnd-segStart) < minSegLength
                smallSegs = smallSegs+1;
            end
            tracksIndxTemp = zeros(size(tracksIndxOld));
            tracksIndxTemp(segStart:segEnd-1) = tracksIndxOld(segStart:segEnd-1);
            tracksIndxNew = [tracksIndxNew; tracksIndxTemp];
            tracksCoordTemp = nan(size(tracksCoordOld));
            tracksCoordTemp(((segStart-1)*8+1):(segEnd*8)) = tracksCoordOld(((segStart-1)*8+1):(segEnd*8));
            tracksCoordNew = [tracksCoordNew; tracksCoordTemp];
            totalSegs = totalSegs+1;
            tracksDivided(i).originSegment(iTrack) = j;
            segStart = segEnd;
            iTrack = iTrack+1;
        end
    end
    tracksDivided(i).tracksFeatIndxCG = tracksIndxNew;
    tracksDivided(i).tracksCoordAmpCG = tracksCoordNew;
    tracksDivided(i).seqOfEvents = seqOfEventsNew;
    tracksDivided(i).isInside = isInsideNew;
end

if smallSegs/totalSegs > 0.5
    fprintf(['More than 50%% of segments are smaller than the minimum length \n' ...
    'of %g. Particles may be oscillating between inside and outside at the \n' ...
    'periphery of the masks. Consider increasing the mask size by lowering \n' ...
    'the Gaussian masking threshold or reducing the minimum length argument \n' ...
    'of this function. \n'],minSegLength)
end
    
end

