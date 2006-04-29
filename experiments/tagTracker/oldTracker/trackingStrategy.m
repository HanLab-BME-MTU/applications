function [trackPairs,nsp]=trackingStrategy(idlist);
% TRACKINGSTRATEGY returns frame pairs to be tracked
% 
%
% SYNOPSIS trackPairs=trackingStrategy(idlist);
%
% INPUT idlist         : idlist or idlisttrack
%       options        : strategy options
%
% OUTPUT trackPairs : matrix [source target] numberOfTimepointsx2

% c: 26/3/03	dT
%%%%%%%%%%%%%%%%%%%%%%

% Constants
BW = -1;
FW =+1;
expSpts=4;            % max. number of spots per image (expected) (normal 4, simulation 2)

% General Init 
sl_start=1;
sl_end=length(idlist);  % number of timepoints
trackPairs = [0, 0];

% find entries in 'idlist' where 'expSpts' spots are found
for t=1:sl_end
    if ~isempty(idlist(t).linklist)
        nsp(t)=max(idlist(t).linklist(:,2)); 
    else
        nsp(t)=0;
    end;
end;

% find possible starting frames
%%%possibleSource=findextendedlocmax(nsp);

% For current strategy:


%take all frames before or after which there are more spots
di=sign((nsp(2:end)-nsp(1:end-1))); %use sign, because there could be jumps of more than one, if tp's have been rejected
possibleSource=union(find(di==-1),find(di==1)+1);
%find frames that contain the most spots
maxSpotsIdx = possibleSource(find(nsp(possibleSource) == max(nsp(possibleSource))));

%Remove non-valid starting frames
%only take frames with max # of tags within +/- 10 timepoints, make sure that there
%are no empty starting frames (should not happen, actually)
for i=length(possibleSource):-1:1
    if isempty(idlist(possibleSource(i)).linklist)
        possibleSource(i) = [];
    else %check if it is a global max or an isolated enough local max. if not, remove from list
        if ~any(maxSpotsIdx == possibleSource(i)) & ~any((maxSpotsIdx>possibleSource(i)-10 & maxSpotsIdx<possibleSource(i)+10))
            possibleSource(i) = [];
        end
    end
end;

ctidx=0;
% main loop
for direction= [FW, BW]
    
    %direction-specific init
    for  sourceIdx=possibleSource
        
        %otherSources: other, possibly competing sources: possibleSources that
        %are further down the road (if direction > 0: larger, else: smaller indices)
        otherSources = possibleSource(find(direction * possibleSource > direction * sourceIdx));
        
        %assign otherMax only if sourceIdx is not among maxSpotsIdx!
%         if any(maxSpotsIdx == sourceIdx)
%             otherMax = [];
%         else
            otherMax = maxSpotsIdx(find(direction * maxSpotsIdx > direction * sourceIdx));
%         end
        
        targetIdx=sourceIdx+direction; %targetFrame
        %find valid targetIdx (not empty, not )
        while ((targetIdx>=sl_start & targetIdx<=sl_end ) & isempty(idlist(targetIdx).linklist) )
            targetIdx=targetIdx+direction;
        end
        
        %while: targetIdx is within bounds & tracking from more 2 less spots & targetIdx not a target frame yet & ...
        % targetIdx not closer to another source [& target, starting from a nonMax-source not within 10 of max] & {targetIdx not
        % more than 10 timepoints away from source}
        while( (targetIdx<=sl_end & targetIdx>=sl_start) & (nsp(sourceIdx)> nsp(targetIdx)) & ~any(trackPairs(:,2)==targetIdx) & ...
                   ~any(abs(otherMax - targetIdx) < abs(sourceIdx - targetIdx))) %[~any(abs(otherSources - targetIdx) < abs(sourceIdx - targetIdx)) & ~any(abs(otherMax - targetIdx)<=10)]& {abs(sourceIdx - targetIdx)<=10}
            
            %store current trackPair
            trackPairs=[trackPairs; sourceIdx, targetIdx];
            
            %start looking for next targetIdx for the current source
            targetIdx=targetIdx+direction;
            
            %make sure new target could possibly be valid
            while ( (targetIdx>=sl_start & targetIdx<=sl_end ) & isempty(idlist(targetIdx).linklist)) 
                targetIdx=targetIdx+direction;
            end
            
        end;
    end;
end;
%trackPairs have been initiated with 0,0 to allow logical statement ~any(trackPairs(:,2)==targetIdx)
trackPairs = trackPairs(2:end,:); 
