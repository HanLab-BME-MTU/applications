function [idlist] = LG_deleteFrame(idlist,dataProperties,currentTime,goodTimes,recalc)
%LG_deleteFrame is the function that deletes a frame
% currentFrame could also be a list of frames!

% go ahead and delete
[idlist(currentTime).linklist] = deal([]);

% before currentTime: negative, after: positive
beforeTime = goodTimes(goodTimes - currentTime(1) < 0);
if ~isempty(beforeTime)
    beforeTime = beforeTime(end);
end
afterTime = goodTimes(goodTimes - currentTime(end) > 0);
if ~isempty(afterTime)
    afterTime = afterTime(1);
end

% linkup, linkdown
switch isempty(beforeTime) + 2 * isempty(afterTime)
    case 0 % none empty
        
        idlist(beforeTime).linklist(:,7) = ...
            idlist(afterTime).linklist(:,2);
        idlist(afterTime).linklist(:,6) = ...
            idlist(beforeTime).linklist(:,2);
        
    case 1 % no beforeTime
        
        idlist(afterTime).linklist(:,6) = 0;
        
    case 2 % no afterTime
        
        idlist(beforeTime).linklist(:,7) = 0;
        
    case 3 % this was (sadly) the last frame
        
        h = warndlg('You just deleted the last frame!',...
            'idlist is now empty!');
        uiwait(h);
        
end

% update goodTimes here!
% setxor returns what is not in both. Since currentTime could contain
% non-goodTimes frames, we have to remove those first.
goodTimes = setxor(goodTimes,intersect(goodTimes,currentTime(:)));


% make sure there aren't tags that never belong to a spot
idlist = LG_checkDeadTags(idlist,goodTimes);

% if recalc has been selected, do so
if recalc
    idlist = LG_recalcIdlist(idlist,dataProperties,{0});
end

% remember that we have deleted a frame
ans = {'without recalc';'with recalc'};
idlist(1).stats.status{end+1} = ...
    sprintf('%s: deleted frame %i:%i %s',...
    date,currentTime(1),currentTime(end),ans{recalc+1});
        