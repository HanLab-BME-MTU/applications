function [clustHistoryAll,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState(tracksAggregStateDef)
%CLUSTERHISTORYFROMCOMPTRACKS_AGGREGSTATE determines the size and lifetime of all clusters that formed during a simulation
%
%   The function uses the information conatined in seqOfEvents and
%   aggregState, two fields from the output of aggregStateFromCompTracks,
%   in the field defaultFormatTracks.  
%   clusterHistoryFromCompTracks_aggregStateALT, is another implementation
%   of this function which uses aggregState from the field
%   alternativeFormatTracks but seqOfEvents from defaultFormatTracks.     
%
%   SYNOPSIS: [clustHistoryAll,clustHistoryMerged] = ...
%    clusterHistoryFromCompTracks_aggregState(tracksAggregStateDef)
%
%   INPUT:
%       tracksAggregStateDef:  
%                         the structure of track information including 
%                         aggregState in default format.
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
%   Robel Yirdaw, 09/19/13
%       modified, 02/20/14
%       modified, 04/08/14               
%       modified, 11/18/14
%       modified, 28 May 2015 (Khuloud Jaqaman)
%

    %Determine number of compTracks generated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %040814 - now passing defaultFormatTracks directly for memory reasons.
    [numCompTracks,~] = size(tracksAggregStateDef);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %Cluster history from all compTracks will be saved
    clustHistoryAll = cell(numCompTracks,1);
    
    %For each compTrack
    for compTrackIter=1:numCompTracks  
        %seqOfEvents for current compTrack
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %040814 - now passing defaultFormatTracks directly for memory reasons.        
        seqOfEvents = tracksAggregStateDef(compTrackIter,1).seqOfEvents;                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %111814
        %For experimental data, the first iteration (frame) in seqOfEvents
        %can be some value other than 1.  This is occurs becuase an object
        %can appear or disappear at any time while in simulation, all of
        %the objects are present at frame 1. This affects how aggregState
        %is accessed since the frame # is supposed to correspond to a
        %column in aggregState. So, need to shift frame numbers in
        %seqOfEvents, if necessary, when accessing aggregState. The very
        %first frame # is the shift amount.
        frameShift = seqOfEvents(1,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %091813
        %aggregState will be used to get cluster sizes 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %040814 - now passing defaultFormatTracks directly for memory reasons.        
        aggregState = tracksAggregStateDef(compTrackIter,1).aggregState;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Iteration points of first and last events
        firstEventIter = find(~isnan(seqOfEvents(:,4)),1,'first');
        lastEventIter = find(~isnan(seqOfEvents(:,4)),1,'last');
        %Number of all events
        totClustChngs = (lastEventIter-firstEventIter+1) + ...
            sum(seqOfEvents(firstEventIter:lastEventIter,2)==1)+1;
        %Initialize
        %021714 - added column #7 for size of resulting cluster
        %021814 - added column #8 for flag for association events
        clustHistoryTemp = NaN(totClustChngs,8);
        clustHistIndx = 1;
        
        %Process events in seqOfEvents
        for eventIndx=firstEventIter:lastEventIter         
            %fprintf('\n%d',eventIndx);
            %If event is a dissociation
            if (seqOfEvents(eventIndx,2) == 1)
                %The dissociating receptor's (cluster's) segment number
                clustHistoryTemp(clustHistIndx,1) = seqOfEvents(eventIndx,3);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %091713    
                %Getting cluster size from entries in aggregState

                %Use the row to access seqOfEvents and obtain
                %the segment number in column 3 which corresponds to the
                %newly split segment.  Then access aggregState with this
                %number as the row and the current iteration value (column
                %1 in seqOfEvents of the same row) as the column to
                %get the size of the cluster that split
                %
                %Modified 111814 to accomodate seqOfEvents that start at
                %frames > 1
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEvents(eventIndx,3),...
                    seqOfEvents(eventIndx,1) - frameShift + 1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Starting iteration point of this cluster
                clustHistoryTemp(clustHistIndx,3) = seqOfEvents(eventIndx,1);
                %Increment saved event count
                clustHistIndx = clustHistIndx + 1;
            else
                %Event is association
                %Find the cluster that is ending in clustHistory.  If it
                %doesn't exist, then this event is at the start of the
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
                    
                    %021714 - added column #7 for size of resulting cluster
                    %Modified 111814 to accomodate seqOfEvents from
                    %experimental data which can have NaN entries for
                    %column #4 corresponding to segments appearing or
                    %disappearing. The segment in column #4 is unknown if
                    %NaN, so in this case, the final size is unknown.
                    if (~isnan(seqOfEvents(eventIndx,4)))
                        clustHistoryTemp(endingTrackIndx,7) = ...
                            aggregState(seqOfEvents(eventIndx,4),...
                            seqOfEvents(eventIndx,1) - frameShift + 1);
                    end
                end

            end %If event 1 or 2

            %The cluster from/with which the above is occuring is also
            %changing. Update accordingly.
            %Find the changing segment in clustHistory.  If it doesn't exist, then
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
                
                %021714 - added column #7 for size of resulting cluster
                %Modified 111814
                clustHistoryTemp(changingTrackIndx,7) = ...
                    aggregState(seqOfEvents(eventIndx,4),...
                    seqOfEvents(eventIndx,1) - frameShift + 1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %02/18/14
                %Special case for cluster size 1
                %If this is a dimerization and both monomers are listed,
                %i.e. the start time for both is known and are now merging,
                %then set flag (col. 8) for the ending track to 1 to
                %indicate that this event is listed twice.
                %{
                if (~isempty(endingTrackIndx) &&...
                        (clustHistoryTemp(endingTrackIndx,2) == 1) &&...
                        (clustHistoryTemp(changingTrackIndx,2) == 1) )
                    clustHistoryTemp(endingTrackIndx,8) = 1;
                end
                %} 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                %022014
                %Do the above applied to all cluster sizes - set col.8 to 1 if
                %endingTrack also exists.  Flag is set for both.  This
                %allows for picking associations where start and end time
                %values for both segments are recorded.
                if ((seqOfEvents(eventIndx,2) == 2) && ~isempty(endingTrackIndx) )
                    clustHistoryTemp(endingTrackIndx,8) = 1;
                    clustHistoryTemp(changingTrackIndx,8) = 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %} should really use one of the above blocks.
            end

            %A dissociation or association has occured involving the cluster
            %with segment number on column 4 of seqOfEvents.  This starts a new
            %cluster - save its segment number, current size and starting
            %iteration point in clustHistory.
            %Modified 111814 to accomodate seqOfEvents from experimental 
            %data which can have NaN entries for column #4 corresponding to
            %segments appearing or disappearing. The segment in column #4 
            %is unknown if NaN. So in this case, a new cluster can not be
            %started.
            if (~isnan(seqOfEvents(eventIndx,4)))                      
                %The receptor's (cluster's) track number
                clustHistoryTemp(clustHistIndx,1) = seqOfEvents(eventIndx,4);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %091713  
                %Getting cluster size from entries in aggregState; in this case
                %the row is in column 4 of seqOfEvents since this an update for
                %the changing track            
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEvents(eventIndx,4),...
                    seqOfEvents(eventIndx,1) - frameShift + 1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             

                %Starting iteration point of this cluster
                clustHistoryTemp(clustHistIndx,3) = seqOfEvents(eventIndx,1);    
                %Increment saved event count
                clustHistIndx = clustHistIndx + 1;
            end
            
            %Reset indx variables
            %changingTrackIndx = []; not needed since set on every iter.
            endingTrackIndx = [];

        end  %For each event in seqOfEvents  
        
        %Remove those events that did not terminate by end of simulation.
        clustHistoryTemp(isnan(clustHistoryTemp(:,4)),:) = [];
        
        %Save current clustHistory in collection
        clustHistoryAll{compTrackIter,1} = clustHistoryTemp;
        
        %Reset vars
        clear clustHistoryTemp tracksCoordAmpCG firstEventIter lastEventIter
        
    end %For each compTrack
    
    %combine cluster histories for all tracks into one output variable
    clustHistoryMerged = cat(1,clustHistoryAll{:,1});
    
end %Function
