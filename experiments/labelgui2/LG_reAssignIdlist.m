function [idlist, success, dataProperties] = LG_reAssignIdlist(idlist,currentTime, goodTimes, pdValues, assignFuture, recalc,dataProperties);
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
%                     colunm 1: indices to good tags in idlist, i.e. of all
%                     that are not estimates of single occurrences or
%                     secondary fusions.
%                     column 2: tag numbers that should be at the positions
%                     determined by column 1
%                     column 3: tag numbers that should be secondary
%                     fusions to the positions determined by column 1
%
%        assignFuture: 1 if not only currentTime, but also subsequent ones
%                      should be changed
%        recalc:      wheter to recalculate idlist or not
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
suggestRecalc = 0;
nTimepoints = length(idlist);

% loop through tags of linklist and reassign idlist
oldLinklist = idlist(currentTime).linklist;
newLinklist = oldLinklist;
% assign no tags, except for the ones that don't appear in this frame - the
% estimates of the single occurences
newLinklist((oldLinklist(:,5)~=3),4) = 0;



% loop through pdValues. Assign primary tags. Then, loop through pdValues
% again and assign secondary tags (if any). All non-5/3 tags have to be
% present in the pdValues, so there shouldn't be a risk of missing a tag.

% pdValues(:,1) are the indices in oldLinklist of the good spots (all
% except secondary fusions and estimates of single occurences)

% loop through primary tags
for i = 1:size(pdValues,1)

    newTag = pdValues(i,2);
    newSpotIdx = pdValues(i,1);

    % if newTag is 0, delete the spot
    if newTag == 0
        newLinklist(newSpotIdx,2) = -1;
    else
        % replace primary tag - we always do that, because secondary
        % fusions are not among the spotIndices

        % assign new tag. Everything else stays the same
        newLinklist(newSpotIdx,4) = newTag;
    end

end % primary tag loop

% loop through secondary tags
for i = 1:size(pdValues,1)

    newTag = pdValues(i,3);

    % if newTag is 0, we just continue
    if newTag == 0
        % do nothing
    else
        % there are three possibilities where to assign a fusion:
        % 1) this spot was a fusion before, so we just put the tag there
        % 2) the former tag position is still available - make it into a
        %    fusion by copying from the primary spot
        % 3) the tag was moved to a fusion spot, and its original position
        %    has been taken over by another tag - assign the fusion tag
        %    in a row where no tag has been assigned yet. In this case, we
        %    should not continue this for future tags

        % Careful: All this is based on the assumption that there can't be
        % a triple fusion.

        newSpotIdx = pdValues(i,1);
        % find row
        spotRow = find(oldLinklist(:,2) == oldLinklist(newSpotIdx,2));

        if length(spotRow) > 1
            % spot used to be a fusion

            % find good row
            newTagRow = spotRow(oldLinkList(spotRow,3) == 3);

            % assignment is simple
            newLinklist(newTagRow, 4) = newTag;


        elseif newLinklist(newTag,4) == 0
            % old position is still available

            % assign row
            newTagRow = newTag;

            % since we know the spot already, we can go ahead and assign
            % the new line. Don't forget to correct cols 3 and 4!
            newLinklist(newTagRow,:) = newLinklist(spotRow,:);
            newLinklist(newTagRow,[3,4]) = [3, newTag];


        else
            % find a 'random' position for assignment
            zeroCol4 = find(~newLinklist(:,4));

            % assign row
            newTagRow = zeroCol4(1);

            % copy spot row and rewrite cols 3 and 4
            newLinklist(newTagRow,:) = newLinklist(spotRow,:);
            newLinklist(newTagRow,[3,4]) = [3, newTag];

            % set assignFuture to 0. Warn if it actually changes the value
            if assignFuture == 1
                ans = questdlg(...
                    'Cannot reassign future with shifted fusions',...
                    'Error','Continue','Cancel');
                if ~strcmp(ans,'Continue')
                    success = 0;
                    return
                end
            else
                assignFuture = 0;
            end

        end

    end % if newTag == 0

end % secondary tag loop

% loop through idlist and redo Q-matrices and trackInits (trackInits and
% amplitudes might be quite off, but that's what happens if you don't
% recalc) - of course, we loop only if assignFuture. Otherwise, we just do
% it for the current frame.
% There is no need to worry about deleted spots - they just disappear. We
% don't need to remove them explicitly, because the number of tags must
% not change, and therefore, no ll4==0 can exist at this point.
% In the loop we first read the order of the old spots, and we assign the
% new list of tags. Then, we update the Q-matrices and the trackInits.

% find number of tags
nTags = size(oldLinklist,1);


% loop. If no future, only reassign current time. (linkup, linkdown need a
% seperate loop)
currentTimeIdx = find(currentTime == goodTimes);
if assignFuture
    endTimeIdx = length(goodTimes);
else
    endTimeIdx = currentTimeIdx;
end
currentTimeIsFirst = currentTimeIdx == 1;


% revert and write Q separately, because in the meantime, we might be
% deleting tags
[idlist, goodTimes] = linkerRevertQmatrices(idlist,goodTimes);

newTagList = newLinklist(:,4);
qNames = {'detectQ_Pix','trackQ_Pix'};



for t = goodTimes(currentTimeIdx : endTimeIdx)'

    if t == currentTime
        % store newLinklist
        idlist(t).linklist = newLinklist;
    end
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
if recalc
    [idlist,dataProperties,success] = ...
        LG_recalcIdlist(idlist,dataProperties);
else
    success = 1;
end

% remember what was done
idlist(1).stats.status{end+1,1} = ...
    sprintf('%s: relinked idlist in frame %i',date,currentTime);