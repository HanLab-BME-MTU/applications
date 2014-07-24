function [idlist, badIdx] = LG_checkDeadTags(idlist,goodTimes)
%LG_checkDeadTags removes all tags that never belong to a spot

% find constantly dead spots - either estimated or secondary fusions
flagList = catStruct(3,'idlist.linklist(:,[3,5])');
badIdx = find(all(flagList(:,1,:) == 1 | flagList(:,2,:) == 3,3));

if ~isempty(badIdx)
    % remove the bad entries
   [idlist] = LG_deleteTag(idlist,badIdx,goodTimes);
end