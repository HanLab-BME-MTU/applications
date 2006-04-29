function idlist = LG_deleteSpot(idlist,currentTime,spotNumber,dataProperties,goodTimes)
%LG_deleteSpot deletes individual spots within labelgui2

% find out whether we delete the spot (and recalc to estimate new position
% and amplitude), or whether we simply remove the tag (if single occurence)

% maybe it is possible to have the single occurence fused - take care of
% that. Look first for all the rows belonging to the spot to delete. Among
% them, find the tags to delete. These don't need to be deleted again as
% spots.
spot2deleteIdx = find(idlist(currentTime).linklist(:,2) == spotNumber);
t2dIdx = find(idlist(currentTime).linklist(spot2deleteIdx,5) == 2);
tag2deleteIdx = spot2deleteIdx(t2dIdx);
spot2deleteIdx(t2dIdx) = [];
% if we delete a tag, indices may change!
spGreaterTagIdx = spot2deleteIdx > tag2deleteIdx;
spot2deleteIdx(spGreaterTagIdx) = spot2deleteIdx(spGreaterTagIdx) -1;

if ~isempty(tag2deleteIdx)
    idlist = LG_deleteTag(idlist, tag2deleteIdx, goodTimes);
    % remember what was done
    idlist(1).stats.status{end+1,1} = ...
        sprintf('%s: deleted tag (while deleting spot %i) %s',...
        date,spotNumber,sprintf('%i ', tag2deleteIdx));
end



% now that the tag is deleted, remove spots
if ~isempty(spot2deleteIdx)
    idlist(currentTime).linklist(spot2deleteIdx,3) = 1;

    % make sure there aren't tags that never belong to a spot
    [idlist, badIdx] = LG_checkDeadTags(idlist,goodTimes);

    % check what has been deleted
    deletedIdx = (ismember(spot2deleteIdx,badIdx));

    % if nothing is left, we just stop. If everything is left, we just
    % continue. If some is left, we need to adjust indices
    if all(deletedIdx)
        % done. We don't even need to re-estimate
    else
        if all(~deletedIdx)
            % no need for index-adjustmet
        else
            % so, theoretically, there could be multiple deletedIdx, and
            % multiple spotIndices. In case that would occur, the code will
            % crash, and then this comment hopefully helps to pinpoint the
            % cause. ;)
            deletedIdx = find(deletedIdx);
            spGreaterDelIdx = spot2deleteIdx > deletedIdx;
            spot2deleteIdx(spGreaterDelIdx) = ...
                spot2deleteIdx(spGreaterDelIdx) - 1;
            spot2deleteIdx(deletedIdx) = [];
        end

        % remove spotNumber
        oldSpotNumber = idlist(currentTime).linklist(spot2deleteIdx,2);
        idlist(currentTime).linklist(spot2deleteIdx,2) = 0;

        % also remove entries in Q-matrix, change spot numbers of tags and
        % of linkup/linkdown

        % don't remove Q-entries - this will happen automatically
        % idlist = linkerRemoveQEntries(idlist,spot2deleteIdx,currentTime);

        % update spot numbers
        biggerIdx = idlist(currentTime).linklist(:,2)>oldSpotNumber;
        idlist(currentTime).linklist(biggerIdx,2) = ...
            idlist(currentTime).linklist(biggerIdx,2) - 1;

        % update linkup, linkdown
        previousTime = max(goodTimes(goodTimes < currentTime));
        nextTime = min(goodTimes(goodTimes > currentTime));

        if ~isempty(previousTime)
            biggerIdx = idlist(previousTime).linklist(:,7)>oldSpotNumber;
            idlist(previousTime).linklist(biggerIdx,7) = ...
                idlist(previousTime).linklist(biggerIdx,7) - 1;
        end

        if ~isempty(nextTime)
            biggerIdx = idlist(nextTime).linklist(:,6)>oldSpotNumber;
            idlist(nextTime).linklist(biggerIdx,6) = ...
                idlist(nextTime).linklist(biggerIdx,6) - 1;
        end


        % re-estimate positions
        idlist = LG_recalcIdlist(idlist,dataProperties,{2});


        % remember deleting spot only if spot has been deleted
        idlist(1).stats.status{end+1,1} = ...
            sprintf('%s: deleted spot %i in frame %i',...
            date,spotNumber,currentTime);
    end
end

