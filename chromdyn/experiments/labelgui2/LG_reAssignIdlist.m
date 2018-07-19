function [idlist, success, dataProperties] = LG_reAssignIdlist(idlist,currentTime, goodTimes, pdValues, assignFuture, recalc,dataProperties)
%LG_reAssignIdlist reassigns the tags in the idlist
%
% SYNOPSIS    [idlist, success, dataProperties] = ...
%                   LG_reAssignIdlist(idlist,currentTime, goodTimes, ...
%                        pdValues, assignFuture, recalc,dataProperties);
%
% INPUT  idlist:      idlist
%        currentTime: frame where assignment is changed
%        goodTimes:   list of non-deleted frames
%        pdValues:    matrix indicating what should be reassigned
%                     colunm 1: row-indices to good spots in idlist, i.e.
%                     of all that are not estimates or secondary fusions.
%                     column 2: tag numbers that should be at the positions
%                     determined by column 1
%                     column 3: tag numbers that should be secondary
%                     fusions to the positions determined by column 1
%
%        assignFuture: 1 if not only currentTime, but also subsequent ones
%                      should be changed
%        recalc:      wheter to recalculate idlist or not. If 2, estimate.
%        dataProperties: dataProperties
%
% OUTPUT  idlist:     changed idlist
%         success:    1 if everything went alright
%         dataProperties: if not empty, dataProperties have changed (due to
%                         changed expected number of spots)
%
% c: jonas 11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

success = 0;
% suggestRecalc = 0; 
% nTimepoints = length(idlist);


% first, loop through the good spots and assign tags. Set to estimated spot
% if necessary. Also, mark if primary fusion
% Then, loop through all tags. If they haven't been assigned yet, try to
% put them into the row they were before, otherwise search for where they
% were exchanged.
% Finally, re-estimate positions.



% remember old linklist
oldLinklist = idlist(currentTime).linklist;
newLinklist = oldLinklist;
% assign no tags, except for the ones that don't appear in this frame - the
% estimates of the single occurences
newLinklist(oldLinklist(:,5)~=3,4) = 0;



% loop through primary spots
for i = 1:size(pdValues,1)

    newTag = pdValues(i,2);
    newSpotIdx = pdValues(i,1);

    % if newTag is 0, mark spot as estimate. We'll assign a tag to it later
    if newTag == 0
        newLinklist(newSpotIdx,2) = 0;
        newLinklist(newSpotIdx,3) = 1;
    else
        % replace primary tag - we always do that, because secondary
        % fusions are not among the spotIndices

        % assign new tag.
        newLinklist(newSpotIdx,4) = newTag;
        
        % check for fusion/tracked tag
        if pdValues(i,3)
            if any(newLinklist(newSpotIdx,3) == [2,5])
            newLinklist(newSpotIdx,3) = 5;
            else
                newLinklist(newSpotIdx,3) = 3;
            end
        else
            if any(newLinklist(newSpotIdx,3) == [2,5])
                newLinklist(newSpotIdx,3) = 2;
            else
                newLinklist(newSpotIdx,3) = 0;
            end            
        end
    end

end % primary tag loop


% in case a spot disappeared: set correct spot numbers
spotIdx = newLinklist(:,2) > 0;
[dummy, newSpotNumbers] = sort(newLinklist(spotIdx,2));
newLinklist(spotIdx,2) = newSpotNumbers;


% loop through remaining tags (in fact, loop through all) and try to assign
% them to a row

for iTag = 1:size(oldLinklist,1)
    
    % skip if tag has already been assigned to good spot
    if any(newLinklist(:,4) == iTag)
        % skip
    else
        % find a row to place the tag. With a strictly limited number of
        % possible tags, there will always be some kind of circular
        % rearrangement between tags. Therefore, go and look in original
        % tag row. If it's occupied, check the row the occupying tag came
        % from. Continue until an empty row is found
        newTagRow = 0;
        possibleTagRow = iTag;
        while ~newTagRow
            if newLinklist(possibleTagRow,4) == 0
                newTagRow = possibleTagRow;
            else
                possibleTagRow = newLinklist(possibleTagRow,4);
            end
        end
        
        % assing tag
        newLinklist(newTagRow,4) = iTag;
        
        % check for fusion
        fusionIdx = pdValues(:,3) == iTag;
        if any(fusionIdx)
            % read position, spotNumber, chi2 (not amplitude!)
            newLinklist(newTagRow,[2,9:12]) = ...
                newLinklist(pdValues(fusionIdx,1),[2,9:12]);
            % make secondary fusion
            newLinklist(newTagRow,3) = 4;
        else
            % just make estimated spot
            newLinklist(newTagRow,[2,3]) = [0,1];
        end
        
    end % if ~skip
end % loop remaining tags

% debug
if any(newLinklist,4) == 0
    error('exception not handled in LG_reAssignIdlist')
end


% assing linklist
idlist(currentTime).linklist = newLinklist;




% loop through idlist and redo Q-matrices and trackInits (trackInits and
% amplitudes might be quite off, but that's what happens if you don't
% recalc) - of course, we loop only if assignFuture. Otherwise, we just do
% it for the current frame.
% There is no need to worry about deleted spots - they just disappear. We
% don't need to remove them explicitly, because the number of tags must
% not change, and therefore, no ll4==0 can exist at this point.
% In the loop we first read the order of the old spots, and we assign the
% new list of tags. Then, we update the Q-matrices and the trackInits.


% loop. If no future, only reassign current time. (linkup, linkdown need a
% seperate loop)
currentTimeIdx = find(currentTime == goodTimes);
if assignFuture
    endTimeIdx = length(goodTimes);
else
    endTimeIdx = currentTimeIdx;
end
%currentTimeIsFirst = currentTimeIdx == 1;



% revert and write Q separately, because in the meantime, we might be
% deleting tags
[idlist, goodTimes] = linkerRevertQmatrices(idlist,goodTimes);

newTagList = newLinklist(:,4);

for t = goodTimes(currentTimeIdx : endTimeIdx)'

    % add new tagList
    idlist(t).linklist(:,4) = newTagList;

    % update trackInit: for every entry, check whether the spot is still an
    % estimated one, then change the tagIdx
    for i = 1:size(idlist(t).trackInit,1)
        oldTagIdx = idlist(t).trackInit(i,1);
        if idlist(t).linklist(oldTagIdx,3) == 1
            idlist(t).trackInit(i,1) = newTagList(oldTagIdx);
        end
    end
    
     % sort linklist - we don't need any information about the old
     % configuration anymore
    idlist(t).linklist = sortrows(idlist(t).linklist,4);
end

% have a separate loop for linkup, linkdown. It's a bit of an overkill, but
% it's much simpler to code
for t = [goodTimes';NaN,goodTimes(1:end-1)']

    % linkup to t-1, linkdown from t-1
    if ~isnan(t(2))
        % linkdown
        idlist(t(2)).linklist(:,7) = idlist(t(1)).linklist(:,2);
        % linkup
        idlist(t(1)).linklist(:,6) = idlist(t(2)).linklist(:,2);
    end

end

% as it turns out, we're not removing tags within this function, after all.
% write Q-matrices again
idlist = linkerWriteQmatrices(idlist,goodTimes);

% make sure there aren't tags that never belong to a spot
idlist = LG_checkDeadTags(idlist,goodTimes);

% make sure that whatever is connected to a good spot stays a good spot
flagList = catStruct(2,'idlist.linklist(:,5)');
% those we need to change have sometimes 0/1 and sometimes 2/3 as tagFlag
badIdx = find(...
    any(ismember(flagList,[0,1]),2) & any(ismember(flagList,[2,3]),2));
for t = goodTimes'
    changeIdx = idlist(t).linklist(badIdx,5) > 1;
    idlist(t).linklist(badIdx(changeIdx),5) = 0;
end

% if selected, recalc idlist
switch recalc
    case 1
    [idlist,dataProperties,success] = ...
        LG_recalcIdlist(idlist,dataProperties);
    case 2
        [idlist,dataProperties,success] = ...
        LG_recalcIdlist(idlist,dataProperties,2);
    otherwise
    success = 1;
end

% remember what was done
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: relinked idlist in frame %i',date,currentTime);