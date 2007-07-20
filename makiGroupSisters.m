function dataStruct = makiGroupSisters(dataStruct,verbose)
%MAKIGROUPSISTERS groups sister kinetochores
%
% SYNOPSIS: dataStruct = makiGroupSisters(dataStruct)
%
% INPUT dataStruct: data structure as created by makiMakeDataStruct
%       verbose (opt) : 0 - no output (default)
%                       1 - plot 4 frames with sister assignment
%                       2 - 1 & plot all tracks
%
% OUTPUT dataStruct: Same as input, but with field sisterList
%              sisterList is a structure with length equal to the number of
%              sister kinetochore pairs.
%              sisterList(1).trackPairs is an nPairs-by-6 array with
%                   [track1,track2,cost,avg.dist,variance,alignment], that
%                   is sorted according to increasing cost.
%                   track1,2: track indices as in dataStruct.tracks.
%                   cost: cost of grouping
%                   avg.dist: average distance between the two tracks
%                   variance: variance of the distance between the tracks
%                   alignment: f(tan(alpha)), where alpha is the average
%                       angle between the distanceVector and the first
%                       eigenVector of planeFit.eigenVectors
%              sisterList(iPair).distanceVectors is a nTimepoints-by-3 
%                   array with the distance vectors between the tracks (use
%                   normList to get distances and normed vectors). Wherever
%                   the two tracks don't overlap or contain gaps,
%                   distanceVectors contains NaNs
%              sisterList(iPair).coords is a nTimepoints-by-3 array with
%                   the coordinates of the first of the two tracks.
%                   sisterList(i).coords + sisterList(i).distanceVectors
%                   returns the coordinates of the other track.
%
% REMARKS The strategy is to first try and identify sister pairs among the
%           tracks that span at least 75% of the movie
%         Grouping is mainly based on the idea that within sister variance
%           is smaller than between sisters variance. Additionally, there
%           is a user-set distance cutoff above which pairs won't be
%           linked. In metaphase, the alignment of sisters with respect to
%           the plane is taken into account.
%         At the end of the code is a plotting function for the distance
%           between the tracks for debugging
%         The code cannot handle merged/splitted tracks!
%
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 16-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%====================================
%% TEST INPUT & READ PARAMETERS
%====================================

if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
if nargin < 2 || isempty(verbose)
    verbose = 0;
end

% read parameters. 
goodTrackRatio = dataStruct.dataProperties.groupSisters.goodTrackRatio;
maxDist = dataStruct.dataProperties.groupSisters.maxDist;
robust = dataStruct.dataProperties.groupSisters.robust;
costFunction = lower(dataStruct.dataProperties.groupSisters.costFunction);

% read movieLength
nTimepoints = dataStruct.dataProperties.movieSize(4);
% read track statistics. This will work only if no merge/split, i.e. if
% there are only two events per track: a start and a finish
try
    trackStats = catStruct(3,'dataStruct.tracks.seqOfEvents');
catch
    error('makiGroupSisters cannot handle merging/splitting')
end

% check for relative length of tracks
relTrackLength = ...
    squeeze(trackStats(2,1,:)-trackStats(1,1,:)+1)./nTimepoints;
% select above 0.75
goodTracks = find(relTrackLength>=goodTrackRatio);
% count
nGoodTracks = length(goodTracks);



%============================
%% READ TRACK INFORMATION
%============================

% preassign distance/variance matrix. Set zeros here for sparse later
[variances,distances,alignment] = deal(zeros(nGoodTracks));

% read normals to plane - takes time, thus outside the loop
normals = catStruct(2,'dataStruct.planeFit.eigenVectors(:,1)')';

% loop through the good tracks, calculate for every pair median distance
% and variance
if verbose == 2
    figure
    hold on
end
for jTrack = 1:nGoodTracks % loop cols

    % read index of track
    jIdx = goodTracks(jTrack);

    % read track coordinates etc.
    [colCoords,colTime,colIdx] = trackData(jIdx,dataStruct,trackStats);

    % plot individual tracks
    if verbose == 2
        plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),...
            'Color',extendedColors(jTrack))
    end

    for iTrack = jTrack+1:nGoodTracks % loop rows

        % read index of track
        iIdx = goodTracks(iTrack);

        % read track coordinates
        [rowCoords,rowTime,rowIdx] = trackData(iIdx,dataStruct,trackStats);

        % find common time
        [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

        % calculate distance (microns)
        distanceVector = colCoords(colIdx(ctColIdx),:) -...
            rowCoords(rowIdx(ctRowIdx),:);

        % calculate alignment:

        % retain commonTime-normals
        commonNormals = normals(commonTime,:);
        % cos(alpha) = dot(a,b)/(norm(a)*norm(b))
        [distance,distanceVectorN] = normList(distanceVector);
        alpha = acos(abs(dot(distanceVectorN',commonNormals')));
        % average alpha, rather than tan to be nice to pairs that will
        % align eventually. Potentially this can be put under the control
        % of the "robust" switch, too
        meanAlpha = mean(alpha);

        % get average distance, variance
        if robust
            [rMean,rStd]=robustMean(distance);
        else
            rMean=mean(distance);
            rStd = std(distance);
        end
        distances(iTrack,jTrack) = rMean;
        variances(iTrack,jTrack) = rStd^2;
        distances(jTrack,iTrack) = rMean;
        variances(jTrack,iTrack) = rStd^2;
        % alignment: one up to about 30°, then increasing
        %[alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
        %    = deal(tan(meanAlpha)+1-meanAlpha);
        % start at one, increase from the beginning (2 at 45°)
        %[alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
        %    = deal(tan(meanAlpha)+1);
        % start at one, increase from the beginning (2 at 45°)
        %[alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
        %    = deal(tan(meanAlpha)+1);
        % start at one, increase from the beginning (2 at 30°)
        [alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
            = deal(2*sqrt(3)*tan(meanAlpha)+1);

    end
end

%=================================
%% CREATE COST MATRIX & GROUP
%=================================

[r2c,c2r,costMat,linkedIdx] = ...
    linkTracks(distances,variances,alignment,...
    nGoodTracks,maxDist,costFunction);

% get cutoff for visualization
% cutoff distances (b/c they'll be used for refinement
[dummy,cutoffDistance]=cutFirstHistMode(distances(linkedIdx),verbose>0);
if verbose
    set(gcf,'Name',sprintf('Distance cutoff for %s',dataStruct.projectName))
end

goodPairIdxL = r2c==c2r;
if verbose
    % plot for 4 frames (that's about how many can be properly displayed)
    deltaT = max(1,floor(nTimepoints/4));
    tOffset = max(1,ceil(deltaT/2));
    t=tOffset:deltaT:nTimepoints;
    % highlight polygons
    r2cTmp = r2c;
    r2cTmp(~goodPairIdxL) = -r2cTmp(~goodPairIdxL);
    plotGroupResults(t,r2cTmp,nGoodTracks,...
        goodTracks,dataStruct,distances,cutoffDistance,...
        sprintf('Initial grouping for %s. G/B-Cutoff=distance',...
        dataStruct.projectName))
end

if ~strcmp(costFunction,'anaphase')

    % even in metaphase, sister distances seem to be small and tightly
    % distributed. Thus, cutoff with distance and relink (this will not work in
    % anaphase, of course).
    [r2c,c2r,costMat,linkedIdx] = ...
        linkTracks(distances,variances,alignment,...
        nGoodTracks,cutoffDistance,costFunction);


    % check for good pairs (non-polygons)
    goodPairIdxL = r2c==c2r;
    if verbose
        % highlight polygons
        r2cTmp = r2c;
        r2cTmp(~goodPairIdxL) = -r2cTmp(~goodPairIdxL);
        [dummy,cutoff]=cutFirstHistMode(costMat(linkedIdx),0);
        plotGroupResults(t,r2cTmp,nGoodTracks,...
            goodTracks,dataStruct,costMat,cutoff,...
            sprintf('Refined grouping for %s. G/B-Cutoff=cost',...
            dataStruct.projectName))
    end
end % if ~anaphase



%=================================
%% RESOLVE POLYGONS
%=================================

% identify polygons. Polygons have to be closed, thus, it should not matter
% where we start. Also, since the distance between polygons and the rest is
% hopefully fairly large, we don't care about neighborhood.
polygonIdx = find(~goodPairIdxL);
% remove the not-linked tracks
polygonIdx(isnan(r2c(polygonIdx))) = [];
polyList = [];
while ~isempty(polygonIdx)
    polyList(1) = polygonIdx(1);
    polygonIdx(1) = [];
    done = false;
    while ~done
        % look up the row the last corner links to
        nextCorner = r2c(polyList(end));
        % check whether the new corner has already been used
        if any(nextCorner == polyList)
            % if yes, exit. The polygon is complete
            done = true;
        else
            polyList(end+1) = nextCorner; %#ok<AGROW>
            % remove corner from polygonIdx
            polygonIdx(polygonIdx==nextCorner) = [];
        end
    end % identify polygon


    % within the polygon: find closest distance to identify first pair.
    % Remove it, and check for more pairs. This will potentially result in
    % more pairs than removal of large distances starting from a tetragon
    done = false;
    while ~done
        % read current cost matrix
        currentCost = costMat(polyList,polyList);
        currentCost(currentCost==-1) = inf;
        % find pair with lowest cost
        [v,minIdx] = min(currentCost(:));
        [idx1,idx2] = ind2sub(size(currentCost),minIdx);
        % write pair into r2c
        r2c(polyList(idx1)) = polyList(idx2);
        r2c(polyList(idx2)) = polyList(idx1);
        c2r(polyList(idx1)) = polyList(idx2);
        c2r(polyList(idx2)) = polyList(idx1);

        % check whether there are still tracks to link
        polyList([idx1,idx2]) = [];

        if length(polyList) > 1
            % continue
        else
            % clear polyList, write NaN into r2c, c2r
            r2c(polyList) = NaN;
            c2r(polyList) = NaN;
            polyList = [];
            done = true;
        end
    end % resolve individual polygons

end % resolve all polygons

if verbose
    % plot final version
    r2cTmp = r2c;
    r2cTmp(~goodPairIdxL) = -r2cTmp(~goodPairIdxL);
    plotGroupResults(t,r2cTmp,nGoodTracks,...
        goodTracks,dataStruct,costMat,cutoff,...
        sprintf('Final grouping for %s. G/B-Cutoff=cost',...
        dataStruct.projectName))
end



%=====================
%% ASSIGN OUTPUT
%=====================

% sisterList:
%   .trackPairs
%   .coords
%   .distanceVectors
goodPairIdxL = r2c==c2r;
linkedIdx=sub2ind([nGoodTracks nGoodTracks],find(goodPairIdxL),r2c(goodPairIdxL));
nGoodPairs = sum(goodPairIdxL);
sisterList(1:nGoodPairs/2,1) = ...
    struct('trackPairs',[],'coords',NaN(nTimepoints,3),...
    'distanceVectors',NaN(nTimepoints,3));

% write trackPairs. Store: pair1,pair2,cost,dist,var,alignment
sisterList(1).trackPairs = ...
    [goodTracks(goodPairIdxL),goodTracks(r2c(goodPairIdxL)),...
    costMat(linkedIdx),distances(linkedIdx),variances(linkedIdx),...
    alignment(linkedIdx)];
% remove redundancy
sisterList(1).trackPairs(:,1:2) = sort(sisterList(1).trackPairs(:,1:2),2);
sisterList(1).trackPairs = unique(sisterList(1).trackPairs,'rows');
% sort according to cost
sisterList(1).trackPairs = sortrows(sisterList(1).trackPairs,3);

% loop trackPairs to get coords, distance
for iPair = 1:nGoodPairs/2
    [rowCoords,rowTime,rowIdx] = ...
        trackData(sisterList(1).trackPairs(iPair,1),dataStruct,trackStats);
    [colCoords,colTime,colIdx] = ...
        trackData(sisterList(1).trackPairs(iPair,2),dataStruct,trackStats);

    % find common time
    [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

    % calculate distance (microns)
    sisterList(iPair).distanceVectors(commonTime,:) =...
        colCoords(colIdx(ctColIdx),:) -...
        rowCoords(rowIdx(ctRowIdx),:);
    sisterList(iPair).coords(commonTime,:) = rowCoords(rowIdx(ctRowIdx),:);
end % loop goodPairs

dataStruct.sisterList = sisterList;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS & DEBUG HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% link tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r2c,c2r,costMat,linkedIdx] = linkTracks(distances,variances,alignment,nGoodTracks,maxDist,costFunction)
% cutoff distances
distCutoffIdx = distances>maxDist;
distances(distCutoffIdx) = 0;
variances(distCutoffIdx) = 0;
alignment(distCutoffIdx) = 0;

% make cost matrix
switch costFunction
    case 'prophase'
        % variance should work here, too
        costMat = distances.*variances;
    case 'prometaphase'
        costMat = distances.*variances;
    case 'metaphase'
        costMat = distances.*variances.*alignment;
    case 'anaphase'
        costMat = alignment;
end
% replace zeros with -1
costMat(~costMat) = -1;

% lap costMat
[r2c,c2r] = lap(costMat,-1,0,1);

% shorten r2c, c2r. No link is nan
r2c = double(r2c(1:nGoodTracks));
r2c(r2c>nGoodTracks) = NaN;
c2r = double(c2r(1:nGoodTracks));
c2r(c2r>nGoodTracks) = NaN;

linkedIdx=sub2ind([nGoodTracks nGoodTracks],1:nGoodTracks,r2c(1:nGoodTracks)');
linkedIdx(isnan(linkedIdx)) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read track coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coords,time,coordIdx] = trackData(idx,dataStruct,trackStats)

% read track coordinates for j
coords = reshape(dataStruct.tracks(idx).tracksCoordAmpCG,8,[])';
coords = coords(:,1:3);
% read timepoints of the track
time = (trackStats(1,1,idx):trackStats(2,1,idx))';
% remove gaps
time(all(isnan(coords),2)) = [];
% remember indices into colCoords that correspond to the timepoints
coordIdx = time - trackStats(1,1,idx) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT/DEBUG FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotGroupResults(tt,r2c,nGoodTracks,goodTracks,dataStruct,variances,cutoff,figureName)
% plot sister-links for one frame. Number indicates goodTrackIdx
% green: variance below cutoff
% blue : variance above cutoff
% red: track with no partner

if nargin < 7 || isempty(cutoff)
    cutoff = inf;
end
figure('Name',figureName)
ntt = length(tt);
rows = ceil(sqrt(ntt));
cols = ceil(ntt/rows);
for p=1:ntt
    t = tt(p);
    [c2pg,c2pb,c2pr]=deal(nan(3*nGoodTracks,3));

    for i=nGoodTracks:-1:1
        idx1 = goodTracks(i);
        % check out of bounds
        if t<dataStruct.tracks(idx1).seqOfEvents(1,1) || t>dataStruct.tracks(idx1).seqOfEvents(2,1)
            c1 = nan(1,3);
        else
            % read via idx, not time
            tIdx = t-dataStruct.tracks(idx1).seqOfEvents(1,1)+1;
            c1=dataStruct.tracks(idx1).tracksCoordAmpCG((tIdx-1)*8+1:(tIdx-1)*8+3);
        end
        if isnan(r2c(i)) || abs(r2c(i)) > nGoodTracks
            c2=nan(1,3);
            v=NaN;
        else
            idx2 = goodTracks(abs(r2c(i)));
            if t<dataStruct.tracks(idx2).seqOfEvents(1,1) || t>dataStruct.tracks(idx2).seqOfEvents(2,1)
                c2 = nan(1,3);
            else
                % read via idx, not time
                tIdx = t-dataStruct.tracks(idx2).seqOfEvents(1,1)+1;
                c2=dataStruct.tracks(idx2).tracksCoordAmpCG((tIdx-1)*8+1:(tIdx-1)*8+3);
            end
            if r2c(i)<0
                v = nan; % polygons are also red
            else
                v = variances(i,r2c(i));
            end
        end

        if v < cutoff
            c2pg(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        elseif isnan(v)
            c2pr(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        else
            c2pb(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        end


    end
    subplot(rows,cols,p)
    plot3(c2pg(:,1),c2pg(:,2),c2pg(:,3),'*-g')
    hold on
    plot3(c2pb(:,1),c2pb(:,2),c2pb(:,3),'*-b')
    plot3(c2pr(:,1),c2pr(:,2),c2pr(:,3),'*-r')
    text(c2pg(1:3:end,1),c2pg(1:3:end,2),c2pg(1:3:end,3),num2str((1:nGoodTracks)'))
    text(c2pb(1:3:end,1),c2pb(1:3:end,2),c2pb(1:3:end,3),num2str((1:nGoodTracks)'))
    text(c2pr(1:3:end,1),c2pr(1:3:end,2),c2pr(1:3:end,3),num2str((1:nGoodTracks)'))
    grid on
    title(sprintf('Frame %i',t))
end

%% PLOT DISTANCE BETWEEN TRACKS
function distance = plotTrackDistance(iIdxT,jIdxT,goodTracks,dataStruct,trackStats) %#ok<DEFNU>
% plots distance between two tracks vs. time; two connected tracks

jIdx = goodTracks(jIdxT);
iIdx = goodTracks(iIdxT);

[rowCoords,rowTime,rowIdx] = trackData(iIdx,dataStruct,trackStats);
[colCoords,colTime,colIdx] = trackData(jIdx,dataStruct,trackStats);
% find common time
[commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

% calculate distance (microns)
distance = sqrt(sum((colCoords(colIdx(ctColIdx),:) -...
    rowCoords(rowIdx(ctRowIdx),:)).^2,2));

figure('Name',sprintf('tracks %i (r) & %i (b)',iIdxT,jIdxT))

subplot(1,2,1)
plot(commonTime,distance)
ylim([0,max(5,max(distance))])

% also plot 3d track
subplot(1,2,2)
plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),'b')
hold on
plot3(rowCoords(:,1),rowCoords(:,2),rowCoords(:,3),'r')
cc = reshape([colCoords(colIdx(ctColIdx),:),...
    rowCoords(rowIdx(ctRowIdx),:),nan(length(distance),3)]',3,[])';
plot3(cc(:,1),cc(:,2),cc(:,3),'g')
plot3(cc(1:2,1),cc(1:2,2),cc(1:2,3),'*')






