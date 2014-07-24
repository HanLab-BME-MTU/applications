function clustHistoryAll = clusterHistoryFromCompTracks_aggregStateALT(compTracks,compTracksALT)
%CLUSTERHISTORYFROMCOMPTRACKS_AGGREGSTATEALT determines the size and lifetime 
%of all clusters that formed during a simulation.
%
%   The function uses the information conatined in seqOfEvents and
%   aggregState, two fields from the output of aggregStateFromCompTracks,
%   in the field defaultFormatTracks.  Here, events are processed from the
%   default formatted seqOfEvents but cluster sizes from aggregState in
%   alternativeFormatTracks. Processing events using default formatted
%   seqOfEvents as opposed to alternative seqOfEvents preferred because the
%   latter produces a large number of new tracks that will not be found in
%   the primary variables of the simulation recept2clustAssign and
%   clust2receptAssign.
%   clusterHistoryFromCompTracks_aggregState is another implementation
%   of this function which uses aggregState and seqOfEvents from 
%   defaultFormatTracks.     
%
%   INPUT:
%       compTracks:    the structure varaiable returned as one of the fields
%                      of receptorInfoLabeled, the main output of the
%                      simulation program receptorAggregationSimple.
%
%       comTracksALT:  the structure varaiable returned by
%                      aggregStateFromCompTracks.
%
%   OUTPUT:
%       clustHistoryAll:  a 1D cell with rows = number of tracks in
%                         compTracks. Each entry contains a clusterHistory 
%                         table a track in compTracks.  clusterHistory is
%                         a 2D array with each row corresponding to an
%                         association or a dissociation event. The 6 colums
%                         give the following information:
%                         1) Track number
%                         2) Cluster size    
%                         3) Starting iteration
%                         4) Ending iteration
%                         5) Lifetime
%                         6) Event that ended the cluster 
%                           (1 = dissociation, 2 = association)
%
%   Robel Yirdaw, 09/19/13
%

    %Determine number of compTracks generated
    [numCompTracks,~] = size(compTracks);
     
    %Cluster history from all compTracks will be saved
    clustHistoryAll = cell(numCompTracks,1);
    
    %For each compTrack
    for compTrackIter=1:numCompTracks  
        %seqOfEvents and tracksCoordAmpCG for current compTrack
        seqOfEvents = compTracks(compTrackIter,1).seqOfEvents;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %091713
        %seqOfEvents, aggregState and an index correspondence variable
        %from alternative compTracks
        %These will be used to get cluster sizes which are stored in
        %aggregState.
       
        seqOfEventsALT = compTracksALT.alternativeFormatTracks(compTrackIter,1).seqOfEvents;
        aggregState = compTracksALT.alternativeFormatTracks(compTrackIter,1).aggregState;
        alt2defSegCorr = compTracksALT.alternativeFormatTracks(compTrackIter,1).alt2defSegmentCorrespond;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Iteration points of first and last events
        firstEventIter = find(~isnan(seqOfEvents(:,4)),1,'first');
        lastEventIter = find(~isnan(seqOfEvents(:,4)),1,'last');
        %Number of all events
        totClustChngs = (lastEventIter-firstEventIter+1) + ...
            sum(seqOfEvents(firstEventIter:lastEventIter,2)==1)+1;
        %Initialize
        clustHistoryTemp = NaN(totClustChngs,6);
        clustHistIndx = 1;
        
        %Process events in seqOfEvents
        for eventIndx=firstEventIter:lastEventIter         

            %If event is a dissociation
            if (seqOfEvents(eventIndx,2) == 1)
                %The dissociating receptor's (cluster's) track number
                clustHistoryTemp(clustHistIndx,1) = seqOfEvents(eventIndx,3);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %091713    
                %Getting cluster size from entries in aggregState
                %The current event can be at a different row in the
                %alternative seqOfEvents
                row1InSeqOfEveALT = find( (seqOfEventsALT(:,1) == seqOfEvents(eventIndx,1)) &...
                    (seqOfEventsALT(:,2) == 1) & (seqOfEventsALT(:,3) == seqOfEvents(eventIndx,3)) ); 
                %Use the row to access alternative seqOfEvents and obtain
                %the track number in column 3 which corresponds to the
                %newly split track.  Then access aggregState with this
                %number as the row and the current iteration value (column
                %1 in alt. seqOfEvents of the same row) as the column to
                %get the size of the cluster that split
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEventsALT(row1InSeqOfEveALT,3),seqOfEventsALT(row1InSeqOfEveALT,1));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Starting iteration point of this cluster
                clustHistoryTemp(clustHistIndx,3) = seqOfEvents(eventIndx,1);
                %Increment saved event count
                clustHistIndx = clustHistIndx + 1;
            else
                %Event is association
                %Find the cluster that is ending in clustHistory.  If it
                %doesn't exists, then this event is at the start of the
                %simulation and its beginning time is unknown. Only the
                %cluster it is joining will be saved below as a new cluster.
                endingTrackIndx = find(clustHistoryTemp(:,1)==seqOfEvents(eventIndx,3),1,'last');
                if (~isempty(endingTrackIndx))
                    %Ending iteration point of this cluster
                    clustHistoryTemp(endingTrackIndx,4) = seqOfEvents(eventIndx,1);
                    %Cluster lifetime
                    clustHistoryTemp(endingTrackIndx,5) = ...
                        clustHistoryTemp(endingTrackIndx,4) - clustHistoryTemp(endingTrackIndx,3);
                    %Type of event ending cluster (1=dissoc., 2=assoc.)
                    clustHistoryTemp(endingTrackIndx,6) = seqOfEvents(eventIndx,2);
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %091713    
                %Getting cluster size from entries in aggregState (in this
                %case for the changing track treated below)
               
                %The current event can be at a different row in the
                %alternative seqOfEvents
                row1InSeqOfEveALT = find( (seqOfEventsALT(:,1) == seqOfEvents(eventIndx,1)) &...
                    (seqOfEventsALT(:,2) == 2) & (seqOfEventsALT(:,3) == seqOfEvents(eventIndx,3)));
                if (isempty(row1InSeqOfEveALT(:,1)))
                    row1InSeqOfEveALT_Set = find( (seqOfEventsALT(:,1) == seqOfEvents(eventIndx,1)) &...
                        (seqOfEventsALT(:,2) == 2) );
                    if (isempty(row1InSeqOfEveALT_Set))
                        fprintf('\nError in clusterHistoryFromCompTracks_aggregStateALT\n');
                    else
                        row1SearchIter = 1;
                        matchFound = 0;
                        %Search for a match in the indices in seqOfEvents
                        %and alternative seqOfEvents. This is required here
                        %because in the alternative track format, an
                        %association event is recorded as creating a brand
                        %new track, each time.  These new tracks *may* not
                        %be found in the default seqOfEvents, especaially
                        %after a large number of iterations. The index
                        %correspondence variable alt2defSegCorr is used to
                        %do the matching.
                        while (~matchFound && row1SearchIter <= length(row1InSeqOfEveALT_Set))
                            
                            currRow = row1InSeqOfEveALT_Set(row1SearchIter,1);
                            correspondingRow = alt2defSegCorr(:,1) == seqOfEventsALT(currRow,4);
                            if (alt2defSegCorr(correspondingRow,2) == seqOfEvents(eventIndx,4))
                                row1InSeqOfEveALT = currRow;
                                matchFound = 1;
                            else
                                row1SearchIter = row1SearchIter + 2;
                            end
                        end %while 
                    end %If set is empty
                end %If not found initial match
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end %If event 1 or 2

            %The cluster from/with which the above is occuring is also
            %changing. Update accordingly.
            %Find the chaning track in clustHistory.  If it doesn't exist, then
            %this event is at the start of the simulation and its beginning
            %time is unkown.  The newly created cluster will be saved below.
            changingTrackIndx = find(clustHistoryTemp(:,1)==seqOfEvents(eventIndx,4),1,'last');
            if (~isempty(changingTrackIndx))
                %Ending iteration point of this cluster
                clustHistoryTemp(changingTrackIndx,4) = seqOfEvents(eventIndx,1);
                %Cluster lifetime
                clustHistoryTemp(changingTrackIndx,5) = ...
                    clustHistoryTemp(changingTrackIndx,4) - clustHistoryTemp(changingTrackIndx,3);
                %Type of event ending cluster (1=dissoc., 2=assoc.)
                clustHistoryTemp(changingTrackIndx,6) = seqOfEvents(eventIndx,2);
            end

            %A dissociation or association has occured invovling the cluster
            %with track number on column 4 of seqOfEvents.  This starts a new
            %cluster - save its track number, current size and starting
            %iteration point in clustHistory.
            %The receptor's (cluster's) track number
            clustHistoryTemp(clustHistIndx,1) = seqOfEvents(eventIndx,4);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %091713  
            %Getting cluster size from entries in aggregState
            
            %The new cluster that has started as a result of a dissociation
            %or association will have a new track number.  If the event is
            %a dissociation this new track is on column 3 of alt.
            %seqOfEvents, else if an association, it will be on column 4.
            %Get this track number and the current iteration to access
            %aggregState for the size of the current cluster similar to
            %above.

            if (seqOfEvents(eventIndx,2) == 1)
                %If dissociation, we need the second row in alt. sequence
                %of events                
                row2InSeqOfEveALT = row1InSeqOfEveALT + 1;                
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEventsALT(row2InSeqOfEveALT,3),seqOfEventsALT(row2InSeqOfEveALT,1));
            else
                %If association, the first row determined above is used,
                %since both entries share the same track number (on column
                %4 of alt. sequence of events).
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEventsALT(row1InSeqOfEveALT,4),seqOfEventsALT(row1InSeqOfEveALT,1));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
            %Starting iteration point of this cluster
            clustHistoryTemp(clustHistIndx,3) = seqOfEvents(eventIndx,1);    
            %Increment saved event count
            clustHistIndx = clustHistIndx + 1;

        end  %For each event in seqOfEvents  
        
        %Remove those events that did not terminate by end of simulation. 
        clustHistoryTemp(isnan(clustHistoryTemp(:,4)),:) = [];
        
        %Save current clustHistory in collection
        clustHistoryAll{compTrackIter,1} = clustHistoryTemp;
        
        %Reset vars
        clear clustHistoryTemp tracksCoordAmpCG firstEventIter lastEventIter
    end %For each compTrack
    
end %Function
