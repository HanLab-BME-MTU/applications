function [sisterList,trackPairs] = groupSisters(tracks,nTimepoints,poleInfo,varargin)
%MAKIGROUPSISTERS groups sister tracks using maximum weigthed matching
%
% SYNOPSIS: [sisterList,trackPairs] = groupSisters(tracks,nTimerpoints)
%
% INPUT dataStruct: data structure as created by makiMakeDataStruct, with
%                   the fields: "dataProperties", "initCoord", "planeFit" &
%                   "tracks". Field "planeFit" can be empty.
%       nTimepoint: a scalar giving the total number of timepoints in the
%                   movie.
%       poleInfo  : pole information, in the form of movieInfo. Needed only
%                   if useAlignment = 1 or associateSis2Pole = 1.
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
% further modified: Khuloud Jaqaman (October 2012)

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
ip.addParamValue('associateSis2Pole',1,@isscalar);
ip.parse(varargin{:});

% Get optional parameters
verbose=ip.Results.verbose;
maxAngle = ip.Results.maxAngle;
maxDist = ip.Results.maxDist;
minOverlap = ip.Results.minOverlap;
useAlignment = ip.Results.useAlignment;
robust = ip.Results.robust;
associateSis2Pole = ip.Results.associateSis2Pole;

if useAlignment || associateSis2Pole
    assert(~isempty(poleInfo))
end

% select tracks whose length is larger than the minimum overlap
tracksSEL=getTrackSEL(tracks);
goodTracks = find(tracksSEL(:,3)>=minOverlap);
nGoodTracks = length(goodTracks);

%get number of poles in spindle and their coordinates
if useAlignment || associateSis2Pole
    numPole = size(poleInfo(1).xCoord,1);
    poleCoord = zeros(nTimepoints,3,numPole);
    poleCoordTmp = zeros(numPole*nTimepoints,6);
    poleCoordTmp(:,1:2) = vertcat(poleInfo.xCoord);
    poleCoordTmp(:,3:4) = vertcat(poleInfo.yCoord);
    if isfield(poleInfo,'zCoord')
        poleCoordTmp(:,5:6) = vertcat(poleInfo.zCoord);
    end
    poleCoordTmp = poleCoordTmp(:,1:2:end);
    for iPole = 1 : numPole
        poleCoord(:,:,iPole) = poleCoordTmp(iPole:numPole:end,:);
    end
end

%get spindle axis vector
if useAlignment && numPole == 2
    spindleAxisVec = poleCoord(:,:,2) - poleCoord(:,:,1);
    spindleAxisMag = normList(spindleAxisVec);
    spindleAxisVec = spindleAxisVec ./ repmat(spindleAxisMag,1,3);
else
    spindleAxisVec = [];
    useAlignment = 0;
    disp('NOTE: useAlignment disabled because not requested and/or spindle not bipolar')
end

%% READ TRACK INFORMATION

% preassign matrices
[variances,distances,alignment,overlapCost,pairCands] = deal([]);

% read track coordinates etc.
[coords,coordsStd,amp,ampStd,time,idx] = deal(cell(nGoodTracks,1));
for i=1:nGoodTracks
    [coords{i},coordsStd{i},amp{i},ampStd{i},time{i},idx{i}] = getTrackData(tracks(goodTracks(i)));
end

%loop through the good tracks and calculate grouping cost elements
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
        
        % calculate distance
        distanceVector = coords{jTrack}(idx{jTrack}(ctColIdx),:) -...
            coords{iTrack}(idx{iTrack}(ctRowIdx),:);
        [distance,distanceVectorN] = normList(distanceVector);
        
        if useAlignment
            
            %get the angle between distance vector and spindle axis
            numDim = size(spindleAxisVec,2);
            distanceDotAxis = sum(distanceVectorN(:,1:numDim) .* spindleAxisVec(commonTime,:),2);
            alpha = acos(abs(distanceDotAxis));
            
            % average alpha, rather than tan to be nice to pairs that will
            % align eventually. Potentially this can be put under the control
            % of the "robust" switch, too
            %average alpha only over frames where there is a plane (the
            %rest are NaN). If none of the frames have a plane, the average
            %will be NaN.
            %also get the standard deviation of alpha
            meanAlpha = nanmean(alpha);
            stdAlpha = nanstd(alpha); %#ok<NASGU>
            
        end
        
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

if ~any(m)
    sisterList = struct('coords1',[],...
        'coords2',[],'sisterVectors',[],'distances',[]);
    trackPairs=[];
    return;
end

%get the median sister distance in order to recalculate the costs and make
%a new assignment in which, instead of favoring the smallest distance, one
%favors the distances closest to the average distance
sisterDistAve = median(distances(m));

costMat = max(abs(distances-sisterDistAve),eps).*variances.*overlapCost;
if useAlignment
    costMat = costMat.*alignment;
end
m=maxWeightedMatching(nGoodTracks,pairCands,1./costMat);

if ~any(m)
    sisterList = struct('coords1',[],...
        'coords2',[],'sisterVectors',[],'distances',[]);
    trackPairs=[];
    return;
end

%% ASSEMBLE SISTER INFORMATION

nGoodPairs = sum(m);
sisterList(1:nGoodPairs,1) = ...
    struct('coords1',NaN(nTimepoints,6),'coords2',NaN(nTimepoints,6),...
    'sisterVectors',NaN(nTimepoints,6),'distances',NaN(nTimepoints,2),...
    'amp1',NaN(nTimepoints,2),'amp2',NaN(nTimepoints,2),...
    'poleAssign12',NaN(nTimepoints,2),...
    'vecFromPole1',NaN(nTimepoints,6),'vecFromPole2',NaN(nTimepoints,6));

% write trackPairs. Store: pair1,pair2,cost,dist,var,alignment
trackPairs = ...
    [goodTracks(pairCands(m,1)),goodTracks(pairCands(m,2)),...
    costMat(m),distances(m),variances(m),alignment(m)];
trackPairs(isnan(trackPairs(:,6)),6) = 0;

% loop over trackPairs to get sister information
validPairs= find(m);
for i=1:numel(validPairs)
    iPair = validPairs(i);
    
    %get information for first sister
    rowCoords = coords{pairCands(iPair,1)};
    rowCoordsStd = coordsStd{pairCands(iPair,1)};
    rowTime  = time{pairCands(iPair,1)};
    rowIdx = idx{pairCands(iPair,1)};
    rowAmp = amp{pairCands(iPair,1)};
    rowAmpStd = ampStd{pairCands(iPair,1)};
    
    %get information for second sister
    colCoords = coords{pairCands(iPair,2)};
    colCoordsStd = coordsStd{pairCands(iPair,2)};
    colTime  = time{pairCands(iPair,2)};
    colIdx = idx{pairCands(iPair,2)};
    colAmp = amp{pairCands(iPair,2)};
    colAmpStd = ampStd{pairCands(iPair,2)};
    
    %find common time between them
    [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);
    
    %store first sister information
    sisterList(i).coords1(commonTime,:) = ...
        [rowCoords(rowIdx(ctRowIdx),:) rowCoordsStd(rowIdx(ctRowIdx),:)];
    sisterList(i).amp1(commonTime,:) = ...
        [rowAmp(rowIdx(ctRowIdx),:) rowAmpStd(rowIdx(ctRowIdx),:)];
    
    %store second sister information
    sisterList(i).coords2(commonTime,:) = ...
        [colCoords(colIdx(ctColIdx),:) colCoordsStd(colIdx(ctColIdx),:)];
    sisterList(i).amp2(commonTime,:) = ...
        [colAmp(colIdx(ctColIdx),:) colAmpStd(colIdx(ctColIdx),:)];
    
    %calculate the vector connecting the two sisters and its std
    sisterVectors = [colCoords(colIdx(ctColIdx),:) - rowCoords(rowIdx(ctRowIdx),:) ...
        sqrt(colCoordsStd(colIdx(ctColIdx),:).^2 + rowCoordsStd(rowIdx(ctRowIdx),:).^2)];
    sisterList(i).sisterVectors(commonTime,:) = sisterVectors;
    
    %calculate the distance between the two sisters and its std
    sisterDist = sqrt(sum(sisterVectors(:,1:3).^2,2));
    sisterDistStd = sqrt(sum((sisterVectors(:,1:3)./repmat(sisterDist,1,3)).^2 .* ...
        sisterVectors(:,4:6).^2,2));
    sisterList(i).distances(commonTime,:) = [sisterDist sisterDistStd];
    
end % loop goodPairs

%% ASSOCIATE SISTERS WITH POLES

if associateSis2Pole
    
    %assign sisters to poles
    switch numPole
        
        case 1 %monopolar spindle
            %assume monotelic attachment
            
            for iSis = 1 : nGoodPairs
                
                %find distance between kinetochores and pole 1
                distSis1 = nanmean(normList(sisterList(iSis).coords1(:,1:3)-poleCoord(:,:,1)));
                distSis2 = nanmean(normList(sisterList(iSis).coords2(:,1:3)-poleCoord(:,:,1)));
                
                %switch sisters if sister 2 is closer
                if distSis2 < distSis1
                    tmp = sisterList(iSis).coords1;
                    sisterList(iSis).coords1 = sisterList(iSis).coords2;
                    sisterList(iSis).coords2 = tmp;
                    sisterList(iSis).sisterVectors = -sisterList(iSis).sisterVectors;
                    tmp = sisterList(iSis).amp1;
                    sisterList(iSis).amp1 = sisterList(iSis).amp2;
                    sisterList(iSis).amp2 = tmp;
                    trackPairs(iSis,[1 2]) = trackPairs(iSis,[2 1]);
                end
                
                %make assignment
                sisterList(iSis).poleAssign12(:,1) = 1;
                
            end
            
        case 2 %bipolar spindle
            %assume amphitelic attachment
            
            for iSis = 1 : nGoodPairs
                
                %find distance between kinetochores and pole 1
                distSis1 = nanmean(normList(sisterList(iSis).coords1(:,1:3)-poleCoord(:,:,1)));
                distSis2 = nanmean(normList(sisterList(iSis).coords2(:,1:3)-poleCoord(:,:,1)));
                
                %switch sisters if sister 2 is closer
                if distSis2 < distSis1
                    tmp = sisterList(iSis).coords1;
                    sisterList(iSis).coords1 = sisterList(iSis).coords2;
                    sisterList(iSis).coords2 = tmp;
                    sisterList(iSis).sisterVectors = -sisterList(iSis).sisterVectors;
                    tmp = sisterList(iSis).amp1;
                    sisterList(iSis).amp1 = sisterList(iSis).amp2;
                    sisterList(iSis).amp2 = tmp;
                    trackPairs(iSis,[1 2]) = trackPairs(iSis,[2 1]);
                end
                
                %make assignment
                sisterList(iSis).poleAssign12(:,1) = 1;
                sisterList(iSis).poleAssign12(:,2) = 2;
                
            end
            
        otherwise %multipolar spindle
            
    end
    
    %calculate vector connecting each kinetochore to its pole
    for iSis = 1 : nGoodPairs
        for iKin = 1 : 2
            jPole = sisterList(iSis).poleAssign12(:,iKin);
            poleCoordTmp = NaN(nTimepoints,3);
            for iPole = 1 : numPole
                iGood = find(jPole==iPole);
                poleCoordTmp(iGood,:) = poleCoord(iGood,:,iPole);
            end
            eval(['vecTmp = sisterList(iSis).coords' num2str(iKin) '(:,1:3) - poleCoordTmp;'])
            eval(['vecTmp = [vecTmp sisterList(iSis).coords' num2str(iKin) '(:,4:6)];'])
            vecTmp(isnan(vecTmp(:,1)),:) = NaN; %#ok<AGROW>
            eval(['sisterList(iSis).vecFromPole' num2str(iKin) ' = vecTmp;'])
        end
    end
    
end

%% SUB-FUNCTION: READ TRACK INFORMATION

function [coords,coordsStd,amp,ampStd,time,coordIdx] = getTrackData(track)

%get start time and end time of track
startTime = track.seqOfEvents(1,1);
endTime = track.seqOfEvents(2,1);
lifeTime = endTime - startTime + 1; %#ok<NASGU>

%get coordinates and their standard deviation
coords = [track.tracksCoordAmpCG(1:8:end)' track.tracksCoordAmpCG(2:8:end)'  track.tracksCoordAmpCG(3:8:end)'];
coordsStd = [track.tracksCoordAmpCG(5:8:end)' track.tracksCoordAmpCG(6:8:end)'  track.tracksCoordAmpCG(7:8:end)'];

%get amplitudes and their standard deviation
amp = track.tracksCoordAmpCG(4:8:end)';
ampStd = track.tracksCoordAmpCG(8:8:end)';

% remove gaps
time = startTime:endTime;
time(all(isnan(coords),2)) = [];

% remember indices into colCoords that correspond to the timepoints
coordIdx = time - startTime + 1;

