function [clustHistoryAll,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState_new(tracksAggregStateDef)
%CLUSTERHISTORYFROMCOMPTRACKS_AGGREGSTATE determines the size and lifetime of all clusters that formed during a simulation
%
%   The function uses the information contained in seqOfEvents and
%   aggregState, two fields from the output of aggregStateFromCompTracks,
%   in the field defaultFormatTracks.  
%
%   SYNOPSIS: [clustHistoryAll,clustHistoryMerged] = ...
%    clusterHistoryFromCompTracks_aggregState(tracksAggregStateDef,infoTime)
%
%   INPUT:
%       tracksAggregStateDef:  
%                         the structure of track information including 
%                         aggregState in default format.
%       infoTime: Structure with fields:
%           .timeStep       : Time between frames/time points. In units of
%                             interest (e.g. s).
%           .sampleStep     : Sampling time step, in same units as
%                             timeStep. Mostly relevant for
%                             simulated data where simulation time step
%                             might be 0.01 s but sampling time step of
%                             interest is e.g. 0.1 s.
%                             Optional. If not input, then sampleStep =
%                             timeStep.
%           .firstLastTP    : Row vector of first and last time points to
%                             use for calculating rates and densities. In
%                             same units as timeStep. 
%                             If only one value is input, it is taken as
%                             the last time point.
%                             If no value is input, then all time points
%                             are used.
%
%   OUTPUT:
%       clustHistoryAll:  a 1D cell with rows = number of tracks in
%                         compTracks. 
%                         Each entry contains a clusterHistory table for a
%                         track in compTracks. In this version cluster 
%                         history is recorded for all clusters in the
%                         obervation time.
%                         clusterHistory is a 2D array with each row
%                         corresponding to an association or a dissociation
%                         event. The 9 colums give the following information:
%                         1) Track segment number.
%                         2) Cluster size.
%                         3) Start time (same units as input timeStep etc).
%                         4) End time (same units as input timeStep etc).
%                         5) Lifetime (same units as input timeStep etc).
%                         6) Event that started the cluster
%                           (1 = dissociation, 2 = association).
%                         7) Event that ended the cluster 
%                           (1 = dissociation, 2 = association).
%                         8) Resulting cluster size.
%                         9) Association flag - 1 indicates the segment
%                            and its partner are both listed and NaN
%                            indicates only the current segment is listed,
%                            i.e. the partner is not listed.
%       
%       clustHistoryMerged: Same information as in clustHistoryAll but with
%                         all cells merged into one 2D array, i.e.
%                         individual track information is lost.
%
%   Robel Yirdaw, 09/19/13
%       modified, 02/20/14
%       modified, 04/08/14               
%       modified, 11/18/14
%       modified, May/June 2015 (Khuloud Jaqaman) 
%       modified, September/October 2016 (Luciana de Oliveira)
%
% *** NOTE (KJ): 
%       As a result of modifications, columns in output have been
%       reshuffled. Thus output is not backward compatible with Robel's
%       functions, unless the latter are modified accordingly.
%

%% Input

%Nothing to do

%% Cluster history collection

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modification Luciana de Oliveira, 10/2016
    
    %If the seqOfEvents has only one track, I need to account this clust
    %size in the clust history
    
    if size(seqOfEvents,1)==2
        
       %Initialize cluster history
        clustHistoryTemp = NaN(1,9);
        
        % track number
        clustHistoryTemp(1,1)=seqOfEvents(1,3);
        
       %cluster size
        clustHistoryTemp(1,2)= aggregState(seqOfEvents(1,1),...
                seqOfEvents(1,1) - frameShift + 1);
             %start time
                 clustHistoryTemp(1,3)=seqOfEvents(1,1);
                  %end time  
        clustHistoryTemp(1,4)=seqOfEvents(2,1);
         %Cluster lifetime
            clustHistoryTemp(1,5) = ...
                clustHistoryTemp(1,4) - clustHistoryTemp(1,3);    
             
            
            clustHistoryTemp = clustHistoryTemp(:,[1:5 9 6:8]); 
       
%Save current clustHistory in collection

clustHistoryAll{compTrackIter,1} = clustHistoryTemp;

    else
    
    %Iteration points of first and last events
    firstEventIter = find(~isnan(seqOfEvents(:,4)),1,'first');
    lastEventIter = find(~isnan(seqOfEvents(:,4)),1,'last');
    
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MODIFICATION LUCIANA DE OLIVEIRA
    %% Part of function I add in to calculate all cluster sizes for all tracks.
    % Now the function has as output a cluster history with all segments
    % history, I need it for the calculations of bootstraping. In the
    % original function clustHistory have only the history for the clusters
    % that are changing. Now it considers also the segments history for the 
    % segments that starts in the first frame and end in the last frame.
    % The function also considers the segments that don't interact with
    % others.
    
    %Determine which are the segments that begin in the first frame
     trackIndexBegin=find(isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
    
    %Determine which are the segments that are changing
     TrackIndexChanging= find(~isnan(seqOfEvents(:,4)));
       
    %Determine which are the segments that finish in the last frame
     trackIndexEnd=find(isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);   
       

      
    %total number of segments in clustHistory
    
    totSegments = length(trackIndexBegin) + (lastEventIter-firstEventIter+1) + ...
    sum(seqOfEvents(firstEventIter:lastEventIter,2)==1);

    %Initialize cluster history
    %021714 - added column #7 for size of resulting cluster
    %021814 - added column #8 for flag for association events
    %2 June 2015 (KJ) - added column #9 for event that started the cluster
    clustHistoryTemp = NaN(totSegments,9);
    clustHistIndxbegin = 1;
 
    %% cluster history for the segments that start in the first frame

for trackIndexNew=trackIndexBegin(1):trackIndexBegin(end)
    % track number
    clustHistoryTemp(clustHistIndxbegin,1) = trackIndexNew;
    % where it starts
    clustHistoryTemp(clustHistIndxbegin,3) = seqOfEvents(trackIndexNew,1);
     
    %where it finishes  
    segmentFinish=find(seqOfEvents(trackIndexNew,3)==seqOfEvents(TrackIndexChanging,3)|seqOfEvents(TrackIndexChanging,4)&seqOfEvents(TrackIndexChanging,2)==2);
   
   if ~isempty(segmentFinish)
   clustHistoryTemp(clustHistIndxbegin,4)=seqOfEvents(TrackIndexChanging(segmentFinish(1)),1)-1;
   else
     segmentFinishLast=find(seqOfEvents(trackIndexNew,3)==seqOfEvents(trackIndexEnd,3));   
    clustHistoryTemp(clustHistIndxbegin,4)=seqOfEvents(TrackIndexChanging(segmentFinishLast(1)),1);
   end
   
   %Getting cluster size from entries in aggregState
    
               clustHistoryTemp(clustHistIndxbegin,2) = ...
                aggregState(seqOfEvents(trackIndexNew,3),...
                seqOfEvents(trackIndexNew,1) - frameShift + 1);
    
            
   %increase in the clustHistIndx
       clustHistIndxbegin = clustHistIndxbegin + 1;
end  %for track index
     


  %% % Initiate the cluster index for the clusters that are changing
    clustHistIndx=clustHistIndxbegin;
    
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
            
            %Type of event starting cluster (1=dissoc., 2=assoc.)
            clustHistoryTemp(clustHistIndx,9) = seqOfEvents(eventIndx,2);
            
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
            
            %The receptor's (cluster's) segment number
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
                       
                            
            %Type of event starting cluster (1=dissoc., 2=assoc.)
            clustHistoryTemp(clustHistIndx,9) = seqOfEvents(eventIndx,2);
            
            %Increment saved event count
            clustHistIndx = clustHistIndx + 1;
        %find if the segments ending in the last frame are also changing in the
%middle of the comptrack
 
 

       %changingTrackIndx = []; not needed since set on every iter.
        endingTrackIndx = [];
        
    end  %For each event in seqOfEvents
    

    end %track inde    
        
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modification Luciana de Oliveira
%% taking the cluster info for the segments that finish in the last frame

%2016/11/01 LRO: I modificate the size of clustHistTemp and then now is
%necessary to remove rows that are not fill with values

%remove rows with all NaNs
iRowNaN= sum(isnan(clustHistoryTemp),2)==9;
clustHistoryTemp(iRowNaN,:)=[];
 
 % 2016/11/01 LRO: find the end of the segment if it is in the tracks
                % that are in the end of the frames.
                
                segmentsNotEnd=find(isnan(clustHistoryTemp(:,4)));
                
                for segIndex=1:length(segmentsNotEnd)
                    
                segmentEnd=find(clustHistoryTemp(segmentsNotEnd(segIndex),1)==seqOfEvents(trackIndexEnd,3));
                clustHistoryTemp(segmentsNotEnd(segIndex),4)=seqOfEvents(trackIndexEnd(segmentEnd(end)),1);
                
                 %Cluster lifetime
            clustHistoryTemp(segIndex,5) = ...
                clustHistoryTemp(segIndex,4) - clustHistoryTemp(segIndex,3);
                end

%reshuffle columns for more intuitive output

 clustHistoryTemp = clustHistoryTemp(:,[1:5 9 6:8]); 
       
       
%Save current clustHistory in collection

clustHistoryAll{compTrackIter,1} = clustHistoryTemp;



    end

end  % for each index
 
%combine cluster histories for all tracks into one output variable

 clustHistoryMerged = cat(1,clustHistoryAll{:,1});
 end


