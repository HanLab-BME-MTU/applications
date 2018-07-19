function idlist = LG_deleteTag(idlist,deleteIdx,goodTimes)
%deleteTag deletes tags from idlist

% the deleteTag-callback checks whether we should recalc or not

% make sure there is something to delete
if isempty(deleteIdx)
    return
end

% sort deleteIdx, so that we won't run into problems with updating
% tagIndices
deleteIdx = sort(deleteIdx,1,'descend');

% loop through goodTimes.
% Aside from simply removing linklist-rows, we need to:
% - remove entries in Q
% - remove entries in trackInit if applicable
% - make secondary spot into primary if applicable
% - change linkup, linkdown if applicable
% - change maxColor
% - adjust spotIdx if applicable

% there is no need to re-sort linklist: the tagIndices are just removed.

% do Qmatrices first
idlist = linkerRemoveQEntries(idlist,deleteIdx,goodTimes);

for t = [goodTimes';NaN,goodTimes(1:end-1)']

    % save spotNumbers and spotFlags of removed tags
    removedSpotNumbers = idlist(t(1)).linklist(deleteIdx,[2,3]);
    % sort removedSpotNumbers so that we have decreasing spotNumbers
    % (otherwise, updating spotNumbers will go wrong)
    removedSpotNumbers = sortrows(removedSpotNumbers,-1);

    % remove entries in linklist
    idlist(t(1)).linklist(deleteIdx,:) = [];

    % if there were actual spots: check for spotIdx, secondary spots,
    % adjust linkup/linkdown
    % if not: remove trackInit-entry

    for i = 1:length(deleteIdx)
        if removedSpotNumbers(i,1) == 0
            % remove trackInit
            trackRemoveIdx = idlist(t(1)).trackInit(:,1)==deleteIdx(i);
            idlist(t(1)).trackInit(trackRemoveIdx,:) = [];
        else
            % check for secondary fusions. Make primary (don't forget if
            % tracked!)
            fusionIdx = idlist(t(1)).linklist(:,2) == removedSpotNumbers(i,1) &...
                idlist(t(1)).linklist(:,3) == 3;
            idlist(t(1)).linklist(fusionIdx,3) = removedSpotNumbers(i,2);

            % update spotNumbers
            biggerSpotNumberIdx = idlist(t(1)).linklist(:,2) > ...
                removedSpotNumbers(i,1);
            idlist(t(1)).linklist(biggerSpotNumberIdx,2) = ...
                idlist(t(1)).linklist(biggerSpotNumberIdx,2) - 1;
        end

        % in any case, adjust tagIndices
        biggerTagNumberIdx = idlist(t(1)).linklist(:, 4) > ...
            deleteIdx(i,1);
        idlist(t(1)).linklist(biggerTagNumberIdx, 4) = ...
            idlist(t(1)).linklist(biggerTagNumberIdx, 4) - 1;
        % and the trackInits. I have no idea why, but it seems for some
        % frames, there's no trackInit
        if ~isempty(idlist(t(1)).trackInit)
            biggerTrackInitNumberIdx = idlist(t(1)).trackInit(:, 1) > ...
                deleteIdx(i,1);
            idlist(t(1)).trackInit(biggerTrackInitNumberIdx, 1) = ...
                idlist(t(1)).trackInit(biggerTrackInitNumberIdx, 1) - 1;
        end
        % and the spotIndices

    end



    % update linkup/linkdown, even if we might not need to: we cannot know
    % whether something changed in the previous frame!
    if ~isnan(t(2))
        % linkdown
        idlist(t(2)).linklist(:,7) = idlist(t(1)).linklist(:,2);
        % linkup
        idlist(t(1)).linklist(:,6) = idlist(t(2)).linklist(:,2);
    end

end

% adjust maxColor. t(1) is the last entry in goodTimes after the loop
idlist(1).stats.maxColor = idlist(t(1)).linklist(end,4);

% adjust labelcolor
idlist(1).stats.labelcolor(deleteIdx) = [];

% remember what was done
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: deleted tags %s',date,sprintf('%i ', deleteIdx));