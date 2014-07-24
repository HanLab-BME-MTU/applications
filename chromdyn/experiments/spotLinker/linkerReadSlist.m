function [idlist, nTimepoints, nSpots, ampList, goodTimes, goodIdx, pix2mu] = linkerReadSlist(slist, constants)
% linkerReadSlist reads data from slist into idlist and other vars


% data is read into idlist(t).linklist-structure
% [1: time, 2: spot#, 3: spotFlag, 4: tagIdx, 5: tagFlag, 
%  6: linkup, 7: linkdown, 8: amp, 9-11: xyz, 12: sigma0 (detector)]

nTimepoints = length(slist);
pix2mu = constants.pix2mu;
% assign idlist, nSpots, ampList (sum of amplitudes), comList (center of
% mass - weights=1). As temporary fields: distMatAmp, distMatXyz,
% sourceIdxList
idlist(1:nTimepoints) = struct(...
    'linklist',[],'info',[],...
    'distMatAmp',[],'distMatXyz',[],'distMat',[],'stats',[],...
    'sourceIdxList',[], 'trackInit',[], 'centroid', []);


[nSpots, ampList, goodIdx] = deal(zeros(nTimepoints,1));

% read slist
% also build ampList, comList
for t=1:nTimepoints
    if ~isempty (slist(t).sp)
        % collect amplitudes, # of spots, coords in microns
        amp = cat(1,slist(t).sp.amp);
        nSpots(t)=length(amp);
        xyz = cat(1,slist(t).sp.cord).*repmat(pix2mu,nSpots(t),1);

        % zero-Column
        zeroCol = zeros(nSpots(t),1);

        % build linklist - don't enter statistics yet
        idlist(t).linklist = [...
            t*ones(nSpots(t),1),...
            [1:nSpots(t)]',...
            zeroCol,...
            zeroCol,...
            zeroCol,...
            zeroCol,...
            zeroCol,...
            amp,...
            xyz,...
            zeroCol];
        % store statistics (stupid history of fieldnames!)
        idlist(t).info = slist(t).statistics;

        % comList, ampList
        idlist(t).centroid = slist(t).COM.*pix2mu;
        ampList(t) = sum(amp);

        % goodIdx
        goodIdx(t) = 1;
    end
end

% make goodIdx logical, and define goodTimes
goodIdx = logical(goodIdx);
goodTimes = find(goodIdx);
