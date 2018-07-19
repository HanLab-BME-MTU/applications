function idlist = LG_deleteSingleOccurences(idlist, goodTimes)
%LG_deleteSingleOccurences is the function for deleting single occurences

% all we need to do here is find the SO in the first good linklist, and
% send them for deletion

tag2deleteIdx = find(ismember(idlist(goodTimes(1)).linklist(:,5),[2,3]));
% remove bad tags
if ~isempty(tag2deleteIdx)
    idlist = LG_deleteTag(idlist, tag2deleteIdx, goodTimes);
end

% remember what was done
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: removed all single occurences',date);


