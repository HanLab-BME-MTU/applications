function [bg,bgh] = bioGraphPlot_MergeSplitEvents(seqOfEventsALT,nodeIDsIN,startSegIN,endSegIN)
%BIOGRAPHPLOT_MERGESPLITEVENTS plots merge split events using biograph and
%biograph objects from the Bioinformatics Toolbox.
%
%   The plot is generated from sequence of events (seqOfEventsALT) which
%   must be in the alternative formats.  Segments are represented as nodes.
%
%   INPUT:
%           seqOfEventsALT:     sequence of events matrix in the
%                               alternative format as found in alternative
%                               format of compTracks.
%           nodeIDsIN:          a cell array of strings which must be of
%                               the same length as the number of segments
%                               to be plotted. The strings will be used to
%                               label each segment on the plot. This is
%                               optional - default strings of the format
%                               'Node #' will be used when not provided.
%           startSegIN:         starting segment index when plotting a
%                               portion of the seqOfEvents. This is
%                               optional
%           endSegIN:           ending segment index when plotting a
%                               portion of the seqOfEvents. This is
%                               optional.
%
%   OUTPUT:
%           bg:                 the constructed biograph object.
%           bgh:                handle to the plotted biograph objects.
%
%   NOTE:   for large (~ > 100 segments) plotting merging and splitting
%           events will be slow and cumbersome. This is best used for small
%           set of events.
%           When providing an ending segment index, the events for that
%           segment will also be plotted.
%
%   Robel Yirdaw, November, 2013.
%       Modified, 11/14/13.
%

    %Initialize and check input parameters
    inputErr = 0;
    nodeIDs = [];
    startSeg = [];
    endSeg = [];
    if (nargin == 0)
        inputErr = 1;
    end
    if (nargin == 2)
        nodeIDs = nodeIDsIN;
    end
    if (nargin == 4)
        startSeg = startSegIN;
        endSeg = endSegIN;
    end

    %If there is no error with input parameters, continue.    
    if (~inputErr)
        %Determine the largest segment/track number.
        maxTrackNum = max(max(seqOfEventsALT(~isnan(seqOfEventsALT(:,4)),3)),...
            max(seqOfEventsALT(~isnan(seqOfEventsALT(:,4)),4)) );
        %Pull out sequence of events.
        msEvents = seqOfEventsALT(~isnan(seqOfEventsALT(:,4)),:);
        %Initialize connection matrix and number of events.
        connMatrix = zeros(maxTrackNum,maxTrackNum);
        numEvents = length(msEvents(:,1));
        
        %Iterate through each event and set up the connection matrix
        %depending on the type of event, merge or split.
        for eventIter=1:2:numEvents

            eventType = msEvents(eventIter,2);            
            if (eventType == 1)
                %Split event
                seg1 = msEvents(eventIter,4);
                seg2 = msEvents(eventIter,3);
                seg3 = msEvents(eventIter+1,3);

                connMatrix(seg1,seg2) = 1;
                connMatrix(seg1,seg3) = 1;
            elseif (eventType == 2)
                %Merge event
                seg1 = msEvents(eventIter,3);
                seg2 = msEvents(eventIter+1,3);
                seg3 = msEvents(eventIter,4);

                connMatrix(seg1,seg3) = 1;
                connMatrix(seg2,seg3) = 1;
            end
        end %for
        
        %If starting and ending segment indices are given, trim the
        %connection matrix and label strings if provided.
        if (~isempty(startSeg) && ~isempty(endSeg))
            largestRelSeg = max([find(connMatrix(endSeg,:) == 1),endSeg]);
            connMatrix = connMatrix(startSeg:largestRelSeg,startSeg:largestRelSeg);
            if (~isempty(nodeIDsIN))
                nodeIDs = nodeIDsIN(startSeg:largestRelSeg,1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Construct biograph object
        bg = biograph(connMatrix,nodeIDs);
        %Plot
        bgh = view(bg);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end %If no input error
    
end %function


    
    
            
        