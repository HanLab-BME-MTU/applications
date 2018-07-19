function [trackPairs,nsp,qVec] = trackingStrategy(inputStruct, trackAll, numberOfSources)
% TRACKINGSTRATEGY returns frame pairs to be tracked
%
% Current strategy  : Only the frames having a number of spots being the
% maximum observed in the movie are retained as possible source frames. For
% each target frame the distance (in time) to all source frames is
% calculated. The source frames with (1) the lowest distance and (2) lowest
% incertitude (Q matrix) are chosen.
%
% SYNOPSIS [trackPairs,nsp,qVec] = trackingStrategy(idlist)
%
% INPUT  inputStruct: either idlist or structure with fields
%                       - nTimepoints
%                       - qVector
%                       - nSpots
%        trackAll   : (opt) whether or not to track all frames [{0}/1]
%        numberOfSources : (opt) number of sources per target frame [{1}]
%
% OUTPUT trackPairs : matrix [source target]
%        nsp        : number of spots
%        qVec       : Vector with [Qxx, Qyy, Qzz] of the detected positions
%
%c: 11/03 michelle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% TEST INPUT
%============

% default: trackAll = 0, numberOfSources = 1;
if nargin < 2 || isempty(trackAll)
    trackAll = 0;
end
if nargin < 3 || isempty(numberOfSources)
    numberOfSources = 1;
end

%===========


%=============================================
% INITIALIZE AND READ INFORMATION FROM IDLIST
%=============================================

% allow idlist as input for backwards compatibility
if isfield(inputStruct, 'linklist')
    idlist = inputStruct;

    nTimepoints = length(idlist);

    film_end = nTimepoints; %last frame of the movie = number of time points

    nsp = (zeros(nTimepoints,1));
    qVec = zeros(nTimepoints,3);

    %fills the vector nsp with the number of spots appearing on each frame
    % and fills vector of uncertainties
    for t = 1:film_end
        if ~isempty(idlist(t).linklist) %if the current time point was not deleted
            nsp(t) = max(idlist(t).linklist(:,2)); %the number of spots is returned
            % qVec: sum of Q-matrix diagonal
            qVec(t,:) = sum(diag(idlist(t).info.detectQ_Pix));
        end;
    end;

else
    nsp = inputStruct.nSpots;
    nTimepoints = inputStruct.nTimepoints;
    qVec = inputStruct.qVector;
end

maxi = max(nsp); %maximum number of spots found during the whole movie

%==============================================

%==============================
% SELECT TARGETS, FIND SOURCES
%==============================

%the frames with a maximum number of spots are used as sources
% the sourceIdx is also the timepoint of the source
sourceIdx = find(nsp == maxi);

% check if we have enough sources - if trackAll, we have one less!
numSources = length(sourceIdx);
if numSources < numberOfSources
    warning('TRACKINGSTRATEGY:NotEnoughSources','Not enough sources. Using %i instead of %i',numSources,numberOfSources);
    numberOfSources = numSources;
end

% choose target frames
if trackAll
    targetIdx = find(nsp); % track all frames
else
    targetIdx = find(nsp ~= maxi & nsp ~= 0); %the frames with a smaller number of spots are used as targets
end

nTarget = length(targetIdx); %number of target frames

% initialize trackPairs
trackPairs = zeros(nTarget * numberOfSources - trackAll * numSources, 2);
pairCt = 0;

for i = 1:nTarget %for each target frame
    diff = abs(targetIdx(i) - sourceIdx); %calculates the distance from the target frame to all the source frames

    %the vector diff is sorted in ascending order (diff_sorted), first by
    %the distance from the target frame, and then the uncertainty, and the
    %indices are returned in diff_index
    [diff_sorted, diff_index] = sortrows([diff,(qVec(sourceIdx))],[1,2]);

    % take only the potential sources! Ohterwise we track all to all.
    diff_index = diff_index(1:numberOfSources);
    
    % prevent source from tracking to itself
    if diff_sorted(1) == 0
        diff_index = diff_index(2:end);
    end

    % Since everything is sorted perfectly fine, we can write sources and
    % targets into the list. There is one fewer source for
    % "source"-targets than there is for "non-source" targets - which has
    % been taken care of when we removed difference 0
    numNewPairs = length(diff_index);

    % with only one source at all, there is nothing that will track the
    % source to itself, so there could be no new pairs
    if numNewPairs > 0
    trackPairs(pairCt + 1:pairCt + numNewPairs, :) = ...
        [sourceIdx(diff_index(1:numNewPairs)),...
        repmat(targetIdx(i),[numNewPairs,1])]; %we will take this one!
    pairCt = pairCt + numNewPairs;
    end

end

