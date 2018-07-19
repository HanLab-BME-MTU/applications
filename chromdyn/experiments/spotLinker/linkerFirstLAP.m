function [idlist, tagIndices, intList, maxTagIdx] = linkerFirstLAP(idlist, nTimepoints, nSpots, t1t2, verbose, intAx)
%linerFirstLAP links all the frames with maxSpots tags


% The maximum distance is 3*sigma. That should be a lot for yeast, but
% sufficiently little for mammalian cells (=> maxCost ~3)

% To ensure continuity, we first link all frames with MAXSPOTS tags (or the
% maximum number of tags). In a second step, all frames with fewer tags are
% linked from several frames with MAXSPOTS tags.

% loop through goodTimes, LAP, and write linkup, linkdown. Reuse t1t2

maxTagIdx = 0;
% tagIndices is going to be used for all frames. Therefore, generate so
% that it will be able to support tags from all frames
[tagIndices] = (zeros(nTimepoints,max(nSpots)));
intList = repmat(NaN,nTimepoints,max(nSpots));

for t = t1t2

    % read nsp
    nsp = nSpots(t);

    % lap without extended testing. Automatically augment, because we could
    % be further than maxDistance!
    [down, up] = lap(idlist(t(1)).distMat,-1,0,1);

    % if we have a number > nsp in the first nsp entries of down or up,
    % there was a birth/death. Therefore, set -1
    down(down>nsp(2)) = -1; % there are maximum nsp(2) points to link to
    up(up>nsp(1)) = -1;
    down = down(1:nsp(1));
    up = up(1:nsp(2));

    idlist(t(1)).linklist(:,7) = down(1:nsp(1));
    idlist(t(2)).linklist(:,6) = up(1:nsp(2));

    % set tagIdx and sort linklist accordingly.
    % If there is a birth, add a new tagIdx.
    % If we're here for the first round, add all new indices, otherwise
    % just write indices for t2
    if t(1) == t1t2(1)
        idlist(t(1)).linklist(:,4) = [1:nsp(1)]';
        maxTagIdx = nsp(1);
        tagIndices(t(1),1:nsp(1)) = idlist(t(1)).linklist(:,4)';
        intList(t(1),idlist(t(1)).linklist(:,4)') = idlist(t(1)).linklist(:,8)';
    end

    oldUp = up > 0;
    newUp = up < 0;

    % write the old tagIndices
    idlist(t(2)).linklist(oldUp,4) = ...
        idlist(t(1)).linklist(idlist(t(2)).linklist(oldUp,6),4);
    newMaxTagIdx = maxTagIdx + sum(newUp);
    % write the new tagIndices
    idlist(t(2)).linklist(newUp,4) = [maxTagIdx + 1 : newMaxTagIdx]';
    % store maxIdx, tagIndices
    maxTagIdx = newMaxTagIdx;
    tagIndices(t(2),1:nsp(2)) = idlist(t(2)).linklist(:,4)';
    intList(t(2),idlist(t(2)).linklist(:,4)') = idlist(t(2)).linklist(:,8)';

end % for t = t1t2 link maxN

if verbose
    plot(intAx,intList,'-+');
end