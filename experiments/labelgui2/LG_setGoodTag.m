function [idlist,dataProperties, dataPropertiesHasChanged] = LG_setGoodTag(idlist,goodTimes,tagIdx,dataProperties,setGood)
%LG_setGoodTag is the function that turns one-time appearance-tags into
%good tags or vice versa
%
% SYNOPSIS [idlist,dataProperties, dataPropertiesHasChanged] =...
%       LG_setGoodTag(idlist,goodTimes,tagIdx,dataProperties,setGood)
% 
% INPUT    idlist
%          goodTimes  indices of non-deleted linklists in idlist
%          tagIdx     idx of tag to turn good or bad
%          dataProperties
%          setGood    1 if make good tag, 0 if make bad tag
%
% OUTPUT   idlist
%          dataProperties
%          dataPropertiesHasChanged 1 if the maximum number of good tags
%               has changed
%
% 5/06 help added by jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if setGood
    % remove flags
    for t = goodTimes'
        idlist(t).linklist(tagIdx,5) = 0;
    end
else
    % put flags back. Hide if estimated position or secondary fusion
    for t = goodTimes'
        if ismember(idlist(t).linklist(tagIdx,3),[1,3])
            idlist(t).linklist(tagIdx,5) = 3;
        else
            idlist(t).linklist(tagIdx,5) = 2;
        end
    end
end

% update maxSpots if necessary
newMax = ...
    nnz(~(idlist(t).linklist(:,5) == 2 | idlist(t).linklist(:,5) == 3));
% if setGood = 1, we might be increasing maxspots, else decreasing
if newMax*setGood > dataProperties.MAXSPOTS*setGood
    dataProperties.MAXSPOTS = newMax;
    dataPropertiesHasChanged = 1;
else
    dataPropertiesHasChanged = 0;
end


% remember what was done
goodOrBad = {'bad','good'};
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: set tag %i ''%s''',date,tagIdx,goodOrBad{setGood+1});