function [idlist, nTimepoints, nSpots, ampList, goodTimes, goodIdx, goodTimesM, t1t2, tagIndices, maxTagIdx, intList] = linkerReadIdlist(idlist, constants)
%linkerReadIdlist does the same job as linkerReadSlist, except that it takes an idlist as input


% read number of timepoints
nTimepoints = length(idlist);
% maxTagIdx: only works like this if recalc is not 1
maxTagIdx = idlist(1).stats.maxColor;

% allocate arrays
[ampList, nSpots] = deal(zeros(nTimepoints,1));
goodIdx = false(nTimepoints,1);
[tagIndices,intList] = deal([]);

% switch according to idlist(1).stats.recalc
% 0: we want to recalc all: remove everything except true spots, their
%       coordinates and amplitudes (trackerSpots are true spots!).
% 1: don't do firstLAP: retain tag assignment of true spots. Remove
%       everything else
% 2: don't do secondLAP: retain all tags. Remove position/amplitude
%       estimates. Retain fusions as spots, but put zero amplitude
% 3: don't change anything aside from the Q-matrices
% 10, 11, 12: as 0, 1, 2

% goodTimesM : always if not 0 - indices of maxTags
% t1t2: consecutive pairs with maxTags
% tagIndices: all tag indices of spots in a frame (1+)
% maxTagIdx : overall maximum tag index (1+)
% intList : all spot intensities in a frame (1+)

recalc = mod(idlist(1).stats.recalc{1},10);

% find goodIdx
for t = 1:nTimepoints
    if ~isempty(idlist(t).linklist)
        goodIdx(t) = true;
        % count spots
        nSpots(t) = max(idlist(t).linklist(:,2));
        % find unique spotIndices and calcuate sum of amplitudes
        spotIdx = idlist(t).linklist(:,2) > 0 & ...
            idlist(t).linklist(:,3) ~= 3;
        ampList(t) = sum(idlist(t).linklist(spotIdx,8));
    end
end

% make sure there are no frames with zero spots in here
goodIdx(nSpots == 0) = false;
goodTimes = find(goodIdx);

% revert Q-matrices
[idlist,goodTimes] = linkerRevertQmatrices(idlist,goodTimes);

for t = goodTimes'
    
    

    % adjust idlists
    switch recalc
        case 0 % do everything

            % remove all non-named spots and all secondary fusions
            ll = idlist(t).linklist;
            badRows = ll(:,2) == 0 | ll(:,3) == 4;
            ll(badRows,:) = [];
            ll(:,[3:7,12]) = 0;
            idlist(t).linklist = ll;

            % empty trackInit
            idlist(t).trackInit = [];

        case 1 % do everything but firstLAP

            % remove all non-named spots, but retain tag-indices. also
            % keep linkup, linkdown
            ll = idlist(t).linklist;
            badRows = ll(:,2) == 0 | ll(:,3) == 4;
            ll(badRows,:) = [];
            ll(:,[3,5,12]) = 0;
            idlist(t).linklist = ll;

            % empty trackInit
            idlist(t).trackInit = [];


        case 2 % do everything after secondLAP

            % remove all estimates of spots and amplitudes (fusions:
            % remove only amplitudes). Also remove all flags
            ll = idlist(t).linklist;
            estIdx = ll(:,3) == 1;
            ll(estIdx,[8:11]) = 0;
            fusIdx = ll(:,3) == 4;
            ll(fusIdx,8) = 0;
            ll(:,12) = 0;
            idlist(t).linklist = ll;
            

            % empty trackInit
            idlist(t).trackInit = [];

        case 3 % do only fusions

            % don't change anything

        otherwise
            error('unrecognized recalc option')
    end
end

% write additional variables
maxNumSpots = min(max(nSpots),constants.MAXSPOTS);
goodTimesM = 0;
delta = 0;
if length(goodTimes) > 2
    while length(goodTimesM) < 2
        goodTimesM = find(nSpots == maxNumSpots-delta);
        delta = delta + 1;
        % make sure we don't get into an infinite loop!
        % Without this clause, this would occur if MAXSPOTS is smaller than
        % max(nSpots), and MAXSPOTS is never found.
        if delta == maxNumSpots 
            if maxNumSpots < max(nSpots)
                currentN = maxNumSpots;
                while currentN < max(nSpots) && length(goodTimesM) < 2
                    currentN = currentN + 1;
                    goodTimesM = find(nSpots == currentN);
                end
            else
                error('there seem to be no spots found at all!')
            end
        end
            
    end
else
    [dummy,goodTimesM] = max(nSpots);
end
t1t2 = [goodTimesM(1:end-1)';goodTimesM(2:end)'];


% for tagIndices and intList: loop again
if recalc == 1
    [tagIndices] = (zeros(nTimepoints,max(nSpots)));
    intList = repmat(NaN,nTimepoints,max(nSpots));
    % we only have as many tags as there are spots in a good frame
    maxTagIdx = nSpots(goodTimesM(1));
    for t = goodTimesM'
        tagIndices(t,1:maxNumSpots) = ...
            idlist(t).linklist(:,4)';
        intList(t,1:maxNumSpots) = ...
            idlist(t).linklist(:,8)';
    end
end

