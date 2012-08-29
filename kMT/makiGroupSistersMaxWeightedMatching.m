function dataStruct = makiGroupSistersMaxWeightedMatching(dataStruct,verbose)
%MAKIGROUPSISTERS groups sister kinetochores
%
% SYNOPSIS: dataStruct = makiGroupSisters(dataStruct)
%
% INPUT dataStruct: data structure as created by makiMakeDataStruct, with
%                   the fields: "dataProperties", "initCoord", "planeFit" &
%                   "tracks". Field "planeFit" can be empty.
%       verbose (opt) : 0 - no plotting (default)
%                       1 - plot 4 frames with sister assignment
%                       2 - 1 & plot all tracks
%
% OUTPUT dataStruct: Same as input, but with field sisterList
%              sisterList is a structure with length equal to the number of
%              sister kinetochore pairs.
%              sisterList(1).trackPairs is an nPairs-by-6 array with
%                   [track1,track2,cost,avg. dist,variance,alignment], that
%                   is sorted according to increasing cost.
%                   track1,2: track indices as in dataStruct.tracks.
%                   cost: cost of grouping
%                   avg. dist: average distance between the two tracks
%                   variance: variance of the distance between the tracks
%                   alignment: f(tan(alpha)), where alpha is the average
%                       angle between the distanceVector and the first
%                       eigenVector of planeFit.eigenVectors
%              sisterList(iPair).coords1 is a nTimepoints-by-6 array with
%                   the coordinates of the first of the two tracks and its
%                   std.
%              sisterList(iPair).coords2 is a nTimepoints-by-6 array with
%                   the coordinates of the second of the two tracks and its
%                   std.
%              sisterList(iPair).sisterVectors is a nTimepoints-by-6 
%                   array with the vector connecting the two sisters and
%                   its std.
%              sisterList(iPair).distances is a nTimepoints-by-2 array with
%                   the distance between sisters and its std.
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
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: Jonas Dorn, Khuloud Jaqaman
% DATE: 16-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TEST INPUT & READ PARAMETERS

if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
if nargin < 2 || isempty(verbose)
    verbose = 0;
end

% read parameters
minOverlap = dataStruct.dataProperties.groupSisters.minOverlap;
if minOverlap < 10
    minOverlap = 10;
end
maxDist = dataStruct.dataProperties.groupSisters.maxDist;
maxAngle = dataStruct.dataProperties.groupSisters.maxAngle * pi / 180;
robust = dataStruct.dataProperties.groupSisters.robust;
useAlignment = dataStruct.dataProperties.groupSisters.useAlignment;
useAnaphase = dataStruct.dataProperties.groupSisters.useAnaphase;

% read movieLength
nTimepoints = dataStruct.dataProperties.movieSize(4);

% read track statistics. This will work only if no merge/split, i.e. if
% there are only two events per track: a start and a finish
try
    trackStats = catStruct(3,'dataStruct.tracks.seqOfEvents');
catch
    error('makiGroupSisters cannot handle merging/splitting')
end


% get track lengths
trackLength = squeeze(trackStats(2,1,:)-trackStats(1,1,:)+1);

% select tracks whose length is larger than the minimum overlap
goodTracks = find(trackLength>=minOverlap);
nGoodTracks = length(goodTracks);

%% READ TRACK INFORMATION

% preassign matrices
[variances,distances,alignment,overlapCost] = deal(NaN(nGoodTracks));

%find frames that have a plane (to calculate alignment cost)
%if none of the frames have a plane, we cannot use the alignment criterion
framesWiPlane = [];
for t = 1 : nTimepoints
    if ~isempty(dataStruct.planeFit) && ~isempty(dataStruct.planeFit(t).planeVectors)
        framesWiPlane = [framesWiPlane; t];
    end
end
if isempty(framesWiPlane)
    useAlignment = 0;
end

%find anaphase frames
if ~isempty(dataStruct.planeFit)
    framePhase = vertcat(dataStruct.planeFit.phase);
else
    framePhase = repmat('e',nTimepoints,1);
end
anaphaseFrames = find(framePhase == 'a');
if isempty(anaphaseFrames)
    lastFrameNotAna = nTimepoints;
else
    lastFrameNotAna = anaphaseFrames(1) - 1;
end

%if the whole movie is in anaphase, there's no point in looking for
%sisters. Exit with an empty sisterList.
if length(anaphaseFrames) == nTimepoints
    sisterList = struct('trackPairs',[],'coords1',[],...
        'coords2',[],'sisterVectors',[],'distances',[]);
    dataStruct.sisterList = sisterList;
    return
end

% read normals to plane
normals = NaN(nTimepoints,3);
if useAlignment == 1 && ~isempty(dataStruct.planeFit)
    normals(framesWiPlane,:) = catStruct(2,'dataStruct.planeFit.planeVectors(:,1)')';
end

[sisterList,trackPairs] = groupSisters(trackStats,nTimepoints,...
    'verbose',verbose,'maxAngle',maxAngle,'maxDist',maxDist,...
    'minOverlap',minOverlap,'useAlignment',useAlignement,'robust',robust);

%% assign output to dataStruct
sisterList(1).trackPairs=trackPairs
dataStruct.sisterList = sisterList;
