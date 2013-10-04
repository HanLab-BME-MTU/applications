function clustHistoryAll = clusterHistoryFromCompTracks(compTracks)
%CLUSTERHISTORYFROMCOMPTRACKS determines the size and lifetime of all
%clusters that formed during a simulation.
%
%   The function uses the information conatined in seqOfEvents and
%   tracksCoordAmpCG, two fields in compTracks, to construct cluster
%   lifetime and sizes.  Here, cluster sizes are determined from
%   tracksCoordAmpCG.  In the alternative implementations of this function,
%   clusterHistoryFromCompTracks_aggregState and
%   clusterHistoryFromCompTracks_aggregStateALT, the variable aggregState
%   is used to determine cluster sizes.
%
%   Note: to determine cluster size using tracksCoordAmpCG, the intensity
%         values need to be rounded to nearest integer.  In some cases,
%         i.e. when the intensity is near half an integer, the rounding can
%         give incorrect cluster sizes. Therefore, better to use the
%         alternative implementaions mentioned above.
%         
%
%   INPUT:
%       comTracks:  the structure varaiable returned as one of the fields
%                   of receptorInfoLabeled, the main output of the
%                   simulation program receptorAggregationSimple.
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
        %Set the following to use intensity amplitudes for cluster sizes
        tracksCoordAmpCG = compTracks(compTrackIter,1).tracksCoordAmpCG;  
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

                %Getting cluster size from intensity amplitude
                clustHistoryTemp(clustHistIndx,2) = ...
                    round(tracksCoordAmpCG(seqOfEvents(eventIndx,3),((seqOfEvents(eventIndx,1)-1)*8+4)));
                
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
            end %If event 1 or 2

            %The cluster from/with which the above is occuring is also
            %changing. Update accordingly.
            %Find the chaning track in clustHistoryTemp.  If it doesn't exist, then
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

            %Get its size from the intensity amplitude in tracksCoordAmpCG
            clustHistoryTemp(clustHistIndx,2) = ...
                    round(tracksCoordAmpCG(seqOfEvents(eventIndx,4),((seqOfEvents(eventIndx,1)-1)*8+4)));        
                    
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
