function [clustHistoryAll,clustHistoryMerged] = ...
    clusterHistoryFromCompTracksComplete(tracksAggregStateAlt)
%clusterHistoryFromCompTracksComplete determines the size and lifetime of all clusters that formed during a simulation
%
%   SYNOPSIS: [clustHistoryAll,clustHistoryMerged] = ...
%    clusterHistoryFromCompTracksComplete(tracksAggregStateAlt)
%
%   INPUT:
%       tracksAggregStateAlt:  
%                         the structure of track information including 
%                         aggregState in alternative format.
%
%   OUTPUT:
%       clustHistoryAll:  a 1D cell with rows = number of tracks in
%                         compTracks. 
%                         Each entry contains a clusterHistory table for a
%                         track in compTracks. Cluster history is only
%                         recorded for clusters that start and end during
%                         the obervation time.
%                         clusterHistory is a 2D array with each row
%                         corresponding to an association or a dissociation
%                         event. The 6 colums give the following information:
%                         1) Track segment number
%                         2) Cluster size    
%                         3) Starting iteration
%                         4) Ending iteration
%                         5) Lifetime
%                         6) Event that ended the cluster 
%                           (1 = dissociation, 2 = association)
%                         7) Resulting cluster size
%                         8) Association flag - 1 indicates the segment
%                            and its partner are both listed and NaN
%                            indicates only the current segment is listed,
%                            i.e. the partner is not listed. 
%       clustHistoryMerged: Same information as in clustHistoryAll but with
%                         all cells merged into one 2D array, i.e.
%                         individual track information is lost.
%
%   Khuloud Jaqaman, May 2015, starting with
%   clusterHistoryFromCompTracks_aggregState
%

%Determine number of compTracks
numCompTracks = length(tracksAggregStateAlt);

%Cluster history from all compTracks will be saved
clustHistoryAll = cell(numCompTracks,1);

%For each compTrack
for compTrackIter=1:numCompTracks
        
    %seqOfEvents for current compTrack
    seqOfEvents = tracksAggregStateAlt(compTrackIter,1).seqOfEvents;
    
    %get aggregation state matrix
    aggregState = full(max(tracksAggregStateAlt(compTrackIter,1).aggregState,[],2));
    numSeg = size(aggregState,1);
    
    clustHistoryTmp = NaN(numSeg,8);
    
    %go over each segment and extract cluster history
    for iSeg = 1 : numSeg
        
        %determine row with information on segment start
        iRowS = find( seqOfEvents(:,3)==iSeg & seqOfEvents(:,2)==1 ... %due to start or split
            | seqOfEvents(:,4)==iSeg & seqOfEvents(:,2)==2 ); %due to merge
        iRowS = iRowS(1);
        startTime = seqOfEvents(iRowS,1) * seqOfEvents(iRowS,4) / seqOfEvents(iRowS,4);
        
        %determine row with information on segment end
        iRowE = find( seqOfEvents(:,3)==iSeg & seqOfEvents(:,2)==2 ... %due to end or merge
            | seqOfEvents(:,4)==iSeg & seqOfEvents(:,2)==1 ); %due to split
        iRowE = iRowE(1);
        endTime = seqOfEvents(iRowE,1) * seqOfEvents(iRowE,4) / seqOfEvents(iRowE,4);
        
        %determine ending event type
        eventType = seqOfEvents(iRowE,2) * seqOfEvents(iRowE,4) / seqOfEvents(iRowE,4);
        
        clustHistoryTmp(iSeg,1:6) = [iSeg aggregState(iSeg) ...
            startTime endTime endTime-startTime eventType];
        
    end
            
    %Save current clustHistory in collection
    clustHistoryAll{compTrackIter,1} = clustHistoryTmp;
    
end %For each compTrack

%combine cluster histories for all tracks into one output variable
clustHistoryMerged = cat(1,clustHistoryAll{:,1});

