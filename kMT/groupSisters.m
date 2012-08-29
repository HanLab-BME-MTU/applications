function [sisterList,trackPairs] = groupSisters(tracks,nTimepoints,varargin)
%MAKIGROUPSISTERS groups sister tracks using maximum weigthed matching
%
% SYNOPSIS: [sisterList,trackPairs] = groupSisters(tracks,nTimerpoints)
%
% INPUT dataStruct: data structure as created by makiMakeDataStruct, with
%                   the fields: "dataProperties", "initCoord", "planeFit" &
%                   "tracks". Field "planeFit" can be empty.
%       nTimepoint: a scalar giving the total number of timepoints in the
%                   movie
%       verbose (opt) : 0 - no plotting (default)
%                       1 - plot 4 frames with sister assignment
%                       2 - 1 & plot all tracks
%
% OUTPUT sisterList is a structure with length equal to the number of
%        sister kinetochore pairs.
%
%        sisterList(iPair).coords1 is a nTimepoints-by-6 array with
%             the coordinates of the first of the two tracks and its std.
%        sisterList(iPair).coords2 is a nTimepoints-by-6 array with
%             the coordinates of the second of the two tracks and its std.
%        sisterList(iPair).sisterVectors is a nTimepoints-by-6
%             array with the vector connecting the two sisters and its std.
%        sisterList(iPair).distances is a nTimepoints-by-2 array with
%             the distance between sisters and its std.
%
%        trackPairs is an nPairs-by-6 array with
%             [track1,track2,cost,avg. dist,variance,alignment], that
%             is sorted according to increasing cost.
%             track1,2: track indices as in dataStruct.tracks.
%             cost: cost of grouping
%             avg. dist: average distance between the two tracks
%             variance: variance of the distance between the tracks
%             alignment: f(tan(alpha)), where alpha is the average
%                  angle between the distanceVector and the first
%                  eigenVector of planeFit.eigenVectors
%
% REMARKS Sister identification is based on globally minimizing (1) the
%           average distance between sisters, (2) the variance of the distance
%           between sisters, and (3) the alignment of sisters with the normal to
%           the metaphase plate (if relevant).
%         Anaphase frames are not used in sister identification.
%         At the end of the code is a plotting function for the distance
%           between the tracks for debugging
%         The code cannot handle merged/splitted tracks!
%
%
% created by: Jonas Dorn, Khuloud Jaqaman
% last modified: Sebastien Besson (May 2012)

%% TEST INPUT & READ PARAMETERS

% Input check
assert(isvector(tracks) && isstruct(tracks(1)));
ip = inputParser;
ip.addOptional('verbose',0,@isscalar);
ip.addParamValue('maxAngle',30 * pi/180,@isscalar);
ip.addParamValue('maxDist',20,@isscalar);
ip.addParamValue('minOverlap',10,@isscalar);
ip.addParamValue('useAlignment',0,@isscalar);
ip.addParamValue('robust',0,@isscalar);
ip.parse(varargin{:});

% Get optional parameters
verbose=ip.Results.verbose;
useAlignment = ip.Results.useAlignment;
maxAngle = ip.Results.maxAngle;
maxDist = ip.Results.maxDist;
minOverlap = ip.Results.minOverlap;
robust = ip.Results.robust;

% select tracks whose length is larger than the minimum overlap
tracksSEL=getTrackSEL(tracks);
goodTracks = find(tracksSEL(:,3)>=minOverlap);
nGoodTracks = length(goodTracks);

%% READ TRACK INFORMATION

% preassign matrices
[variances,distances,alignment,overlapCost, pairCands] = deal([]);

% loop through the good tracks, calculate for every pair mean distance
% and variance


% read track coordinates etc. (one could get the coordinate stds as
% well, but for now no need here - KJ)
% coordinates are in metaphase plane rotated frame of reference
coords = cell(nGoodTracks,1);
time = cell(nGoodTracks,1);
idx = cell(nGoodTracks,1);
coordsStd = cell(nGoodTracks,1);
for i=1:nGoodTracks
    [coords{i},time{i},idx{i},coordsStd{i}] = getTrackData(tracks(goodTracks(i)));
end

iPair=0;
for jTrack = 1:nGoodTracks % loop cols
    % plot individual tracks
    if verbose == 2
        plot(coords{jTrack}(:,1),coords{jTrack}(:,2),'Color',extendedColors(jTrack))
    end
    
    for iTrack = jTrack+1:nGoodTracks % loop rows 
        
        % find common time
        [commonTime,ctColIdx,ctRowIdx] = intersect(time{jTrack},time{iTrack});
        numOverlapFrames = length(commonTime);
        
        %if the common time between the two tracks is at least minOverlap,
        %calculate parameters (otherwise, they stay as NaN, which assigns
        %them -1 in linkTracks)
        if numOverlapFrames < minOverlap, continue; end
        
        
        % calculate distance (microns)
        distanceVector = coords{jTrack}(idx{jTrack}(ctColIdx),:) -...
            coords{iTrack}(idx{iTrack}(ctRowIdx),:);
        
        %get the angle between distance vector and normal
        distance=sqrt(sum(distanceVector.^2,2));
        distanceVectorN = distanceVector./repmat(distance,1,3);
        % [distance,distanceVectorN] = normList(distanceVector);
        alpha = acos(abs(distanceVectorN(:,1)));
        
        % average alpha, rather than tan to be nice to pairs that will
        % align eventually. Potentially this can be put under the control
        % of the "robust" switch, too
        %average alpha only over frames where there is a plane (the
        %rest are NaN). If none of the frames have a plane, the average
        %will be NaN.
        %also get the standard deviation of alpha
        meanAlpha = nanmean(alpha);
        stdAlpha = nanstd(alpha);
        
        if useAlignment && meanAlpha > maxAngle, continue; end
        
        % get distance mean and standard deviation
        if robust
            [rMean,rStd]=robustMean(distance);
        else
            rMean = mean(distance);
            rStd = std(distance);
        end
        
        if rMean > maxDist, continue; end
        
        % Add tracks pair to the list of candidate sisters
        iPair=iPair+1;
        pairCands(iPair,:)= [iTrack,jTrack];
        
        %assign distance mean for pair
        distances(iPair,1) = rMean;
        
        %assign distance variance for pair
        variances(iPair,1) = rStd^2;
        
        %assign alignment cost for pair if the average angle is less
        %than maxAngle degrees. Otherwise, keep as NaN to prohibit the link
        if useAlignment
            alignment(iPair,1) = 2*sqrt(3)*tan(meanAlpha)+1;
        else
            alignment(iPair,1) = NaN;
        end
        
        %assign overlap cost - the longer the overlap, the lower the
        %cost
        overlapCost(iPair,1) = sqrt( 10 / numOverlapFrames );
        
    end %(for iTrack = jTrack+1:nGoodTracks)
end %(for jTrack = 1:nGoodTracks)

%% CREATE COST MATRIX & GROUP

costMat = distances.*variances.*overlapCost;
if useAlignment
    costMat = costMat.*alignment;
end
m=maxWeightedMatching(nGoodTracks,pairCands,1./costMat);
% [r2c,c2r,costMat,linkedIdx] = ...
%     linkTracks(distances,variances,alignment,overlapCost,...
%     nGoodTracks,maxDist,useAlignment);

if ~any(m)
    sisterList = struct('coords1',[],...
        'coords2',[],'sisterVectors',[],'distances',[]);
    trackPairs=[];
    return;
end

%find pairs that get a unique assignment
% goodPairIdxL = r2c==c2r;

% if verbose
%     % plot for 4 frames (that's about how many can be properly displayed)
%     deltaT = max(1,floor(nTimepoints/4));
%     tOffset = max(1,ceil(deltaT/2));
%     t=tOffset:deltaT:nTimepoints;
%     % highlight polygons
%     r2cTmp = r2c;
%     r2cTmp(~goodPairIdxL) = -r2cTmp(~goodPairIdxL);
%     figure;
%     hold on;
%     plotGroupResults(t,r2cTmp,nGoodTracks,...
%         goodTracks,tracks,distances,maxDist,...
%         sprintf('Initial grouping for G/B-Cutoff=distance'))
% end


%% assemble sister information

% sisterList:
%   .coords1
%   .coords2
%   .sisterVectors
%   .distances

% linkedIdx=sub2ind([nGoodTracks nGoodTracks],pairCands(m,1),pairCands(:,2));
nGoodPairs = sum(m);
sisterList(1:nGoodPairs,1) = ...
    struct('coords1',NaN(nTimepoints,6),...
    'coords2',NaN(nTimepoints,6),'sisterVectors',NaN(nTimepoints,6),...
    'distances',NaN(nTimepoints,2));

% write trackPairs. Store: pair1,pair2,cost,dist,var,alignment
trackPairs = ...
    [goodTracks(pairCands(m,1)),goodTracks(pairCands(m,2)),...
    costMat(m),distances(m),variances(m),alignment(m)];
trackPairs(isnan(trackPairs(:,6)),6) = 0;

% % remove redundancy
% trackPairs(:,1:2) = sort(trackPairs(:,1:2),2);
% trackPairs = unique(trackPairs,'rows');

% sort according to cost
% trackPairs = sortrows(trackPairs,3);

% loop over trackPairs to get their coordinates and distances
validPairs= find(m);
for i=1:numel(validPairs)
    iPair = validPairs(i);
    
    %get information for first sister
    rowCoords = coords{pairCands(iPair,1)};
    rowCoordsStd = coordsStd{pairCands(iPair,1)};
    rowTime  = time{pairCands(iPair,1)};
    rowIdx = idx{pairCands(iPair,1)};
    
    %get information for second sister
    colCoords = coords{pairCands(iPair,2)};
    colCoordsStd = coordsStd{pairCands(iPair,2)};
    colTime  = time{pairCands(iPair,2)};
    colIdx = idx{pairCands(iPair,2)};
    
    %find common time between them
    [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);
    
    %store the coordinates of the first sister
    sisterList(i).coords1(commonTime,:) = ...
        [rowCoords(rowIdx(ctRowIdx),:) rowCoordsStd(rowIdx(ctRowIdx),:)];
    
    %store the coordinates of the second sister
    sisterList(i).coords2(commonTime,:) = [colCoords(colIdx(ctColIdx),:) ...
        colCoordsStd(colIdx(ctColIdx),:)];
    
    %calculate the vector connecting the two sisters and its std (microns)
    sisterVectors = [colCoords(colIdx(ctColIdx),:) - rowCoords(rowIdx(ctRowIdx),:) ...
        sqrt(colCoordsStd(colIdx(ctColIdx),:).^2 + rowCoordsStd(rowIdx(ctRowIdx),:).^2)];
    sisterList(i).sisterVectors(commonTime,:) = sisterVectors;
    
    %calculate the distance between the two sisters and its std (microns)
    sisterDist = sqrt(sum(sisterVectors(:,1:3).^2,2));
    sisterDistStd = sqrt(sum((sisterVectors(:,1:3)./repmat(sisterDist,1,3)).^2 .* ...
        sisterVectors(:,4:6).^2,2));
    sisterList(i).distances(commonTime,:) = [sisterDist sisterDistStd];
    
end % loop goodPairs

%% remove extra large distances from sister pairing

%put all sister distances in one vector
% NB: the original makiGroupSisters was slicing the distances using
% lastFramenotAna. Maybe this can be passed as a param/value pair
sisterDist=arrayfun(@(x) x.distances(:,1),sisterList,'Unif',false);
sisterDist = vertcat(sisterDist{:});

% Detect outlier
outlierIndx = detectOutliers(sisterDist,2.5);
[iTime,iPair]=ind2sub([nTimepoints nGoodPairs],outlierIndx);

%remove all distances larger than this maximum distance
for i = 1 : numel(iTime)
    
    
    %remove those timepoints from the sister information
    sisterList(iPair(i)).coords1(iTime(i),:) = NaN;
    sisterList(iPair(i)).coords2(iTime(i),:) = NaN;
    sisterList(iPair(i)).sisterVectors(iTime(i),:) = NaN;
    sisterList(iPair(i)).distances(iTime(i),:) = NaN;
    
end %(for iPair = 1 : nGoodPairs/2)


%% read track coordinates
function [coords,time,coordIdx,coordsStd] = getTrackData(track)

%get indices of feature making track
% featIndx = track.tracksFeatIndxCG;

%get start time and end time of track
startTime = track.seqOfEvents(1,1);
endTime = track.seqOfEvents(2,1);
lifeTime = endTime - startTime + 1;


coords = [track.tracksCoordAmpCG(1:8:end)' track.tracksCoordAmpCG(2:8:end)'  track.tracksCoordAmpCG(3:8:end)'];
coordsStd = [track.tracksCoordAmpCG(5:8:end)' track.tracksCoordAmpCG(6:8:end)'  track.tracksCoordAmpCG(7:8:end)'];

% remove gaps
time = startTime:endTime;
time(all(isnan(coords),2)) = [];

% remember indices into colCoords that correspond to the timepoints
coordIdx = time - startTime + 1;

