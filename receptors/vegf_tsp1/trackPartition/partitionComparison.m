function out = partitionComparison(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Inputs: string name of condition 1, path to CML1, string name of condition 2,
% path to CML2, ...

% Collect input MovieLists into one cell array
% Each cell of this array contains a cell array of information from the
% movies under each condition
MDInfo = cell(1,nargin/2);
conditions = cell(1,nargin/2);
for i = 1:nargin/2
    conditions{i} = varargin{2*i-1};
    file = load(varargin{2*i});
    if isfield(file,'CML')
        MDInfo{i} = loadCML(file.CML);
    elseif isfield(file,'ML')
        MDInfo{i} = loadML(file.ML);
    else
        error('Input %g is not a valid ML or CML path',i)
    end
end

out = cell(1,nargin/2);
for i = 1:nargin/2
    out{i} = parcellfun_progress(@getInfo,MDInfo);
end


end

function MDInfo = loadML(ML)
nMD = numel(ML.movieDataFile_);
MDInfo = cell(1,nMD);
valid = ones(1,nMD);
for i = 1:nMD
    file = load(ML.movieDataFile_{i});
    if file.MD.getProcessIndex('TrackPartitionProcess') > 0
        process = file.MD.processes_{file.MD.getProcessIndex('TrackPartitionProcess')};
    else
        % If MD does not have TrackPartitionProcess, skip it
        valid(i) = 0;
    end
    
    processOut = load(process.outFilePaths_{1});
    MDInfo{i}.tracks = processOut.tracksPart;
    MDInfo{i}.diffAnalysisRes = processOut.diffAnalysisRes;
    MDInfo{i}.info = processOut.info;
end
% Use only MDs with TrackPartitionProcess
MDInfo = MDInfo{valid};
end

function MDInfo = loadCML(CML)
nML = numel(CML.movieListDirectory_);
MDInfo = {};
for i = 1:nML
    file = load(CML.movieListDirectory_{i});
    trackInfoTemp = loadML(file.ML);
    MDInfo = [MDInfo,trackInfoTemp];
end
end

function out = getInfo(MDInfo)

isInside = arrayfun(@(x) x.isInside,MDInfo.tracks);
isOutside = ~isInside;
tracksInside = MDInfo.tracks(isInside);
tracksOutside = MDInfo.tracks(isOutside);
diffInside = MDInfo.diffAnalysisRes(isInside);
diffOutside = MDInfo.diffAnalysisRes(isOutside);

seqInside = arrayfun(@(x) x.seqOfEvents,tracksInside,'UniformOutput',false);
seqInsideMat = cell2mat(seqInside);
seqOutside = arrayfun(@(x) x.seqOfEvents,tracksOutside,'UniformOutput',false);
seqOutsideMat = cell2mat(seqOutside);

nMergeInside = sum((seqInsideMat(:,2) == 1)&(~isnan(seqInsideMat(:,4)))); % need to stop trackPartition from get rid of short merge/splits 
nSplitInside = sum((seqInsideMat(:,2) == 2)&(~isnan(seqInsideMat(:,4))));
nMergeOutside = sum((seqInsideMat(:,2) == 1)&(~isnan(seqInsideMat(:,4)))); 
nSplitOutside = sum((seqInsideMat(:,2) == 2)&(~isnan(seqInsideMat(:,4))));

% Count all tracks (including those inside compound tracks)
nInside = size(seqInsideMat,1)/2;
nOutside = size(seqOutsideMat,1)/2;

lengthInside = zeros(nInside,1);
for i = 1:nInside
    lengthInside(i) = seqInsideMat(2*i,1)-seqInsideMat(2*i-1,1);
end
totalLengthInside = sum(lengthInside);
meanLengthInside = mean(lengthInside);
longestInside = max(lengthInside);

lengthOutside = zeros(nOutside,1);
for i = 1:nOutside
    lengthOutside(i) = seqOutsideMat(2*i,1)-seqOutsideMat(2*i-1,1);
end
totalLengthOutside = sum(lengthOutside);
meanLengthOutside = mean(lengthOutside);
longestOutside = max(lengthOutside);

diffInsideCell = struct2cell(diffInside)';
diffInsideMat = cell2mat(diffInsideCell(:,1));
diffOutsideCell = struct2cell(diffOutside)';
diffOutsideMat = cell2mat(diffOutsideCell(:,1));

asymmetricInside = sum(diffInsideMat(:,1));
notAsymmetricInside = nInside-asymmetricInside;
immobileFullInside = sum(diffInsideMat(:,2) == 0);
confinedFullInside = sum(diffInsideMat(:,2) == 1);
brownianFullInside = sum(diffInsideMat(:,2) == 2);
directedFullInside = sum(diffInsideMat(:,2) == 3);
immobileOneInside = sum(diffInsideMat(:,3) == 0);
confinedOneInside = sum(diffInsideMat(:,3) == 1);
brownianOneInside = sum(diffInsideMat(:,3) == 2);
directedOneInside = sum(diffInsideMat(:,3) == 3);

asymmetricOutside = sum(diffOutsideMat(:,1));
notAsymmetricOutside = nInside-asymmetricOutside;
immobileFullOutside = sum(diffOutsideMat(:,2) == 0);
confinedFullOutside = sum(diffOutsideMat(:,2) == 1);
brownianFullOutside = sum(diffOutsideMat(:,2) == 2);
directedFullOutside = sum(diffOutsideMat(:,2) == 3);
immobileOneOutside = sum(diffOutsideMat(:,3) == 0);
confinedOneOutside = sum(diffOutsideMat(:,3) == 1);
brownianOneOutside = sum(diffOutsideMat(:,3) == 2);
directedOneOutside = sum(diffOutsideMat(:,3) == 3);

[mssSlopeFullInside,genDiffCoefFullInside,scalingPowerFullInside,...
    normDiffCoefFullInside,mssSlopeOneInside,genDiffCoefOneInside,...
    scalingPowerOneInside,normDiffCoefOneInside,confRadiusInside] ...
    = arrayfun(@(x) [x.fullDim.mssSlope,x.fullDim.genDiffCoef,...
    x.fullDim.scalingPower,x.fullDim.normDiffCoef,x.oneDim.mssSlope,...
    x.oneDim.genDiffCoef,x.oneDim.scalingPower,x.oneDim.normDiffCoef,...
    x.confRadInfo.confRadius],diffInside,'UniformOutput',false);

[mssSlopeFullOutside,genDiffCoefFullOutside,scalingPowerFullOutside,...
    normDiffCoefFullOutside,mssSlopeOneOutside,genDiffCoefOneOutside,...
    scalingPowerOneOutside,normDiffCoefOneOutside,confRadiusOutside] ...
    = arrayfun(@(x) [x.fullDim.mssSlope,x.fullDim.genDiffCoef,...
    x.fullDim.scalingPower,x.fullDim.normDiffCoef,x.oneDim.mssSlope,...
    x.oneDim.genDiffCoef,x.oneDim.scalingPower,x.oneDim.normDiffCoef,...
    x.confRadInfo.confRadius],diffOutside,'UniformOutput',false);

% Get tracks that occur pre or post-inside tracks
trackOrigins = arrayfun(@(x) x.originCompoundTrack,MDInfo.tracks);
insideInd = find(isInside);
pre = [];
post = [];
for i = 1:nInside
    ind = insideInd(i); % index of this track within the entire track struct
    origin = tracks(ind).originCompoundTrack;
    siblings = find(trackOrigins == origin); % find other tracks with same origin
    siblingsPre = siblings(siblings < ind); % find sibling tracks before this track
    siblingsPost = siblings(siblings > ind); % find sibling tracks after this track
    
    pre = [pre,siblingsPre];
    post = [post,siblingsPost];
end

diffPre = MDInfo.diffAnalysisRes([pre]);
diffPost = MDInfo.diffAnalysisRes([post]);

diffPreCell = struct2cell(diffPre)';
diffPreMat = cell2mat(diffPreCell(:,1));
diffPostCell = struct2cell(diffPost)';
diffPostMat = cell2mat(diffPostCell(:,1));

asymmetricPre = sum(diffPreMat(:,1));
notAsymmetricPre = nPre-asymmetricPre;
immobileFullPre = sum(diffPreMat(:,2) == 0);
confinedFullPre = sum(diffPreMat(:,2) == 1);
brownianFullPre = sum(diffPreMat(:,2) == 2);
directedFullPre = sum(diffPreMat(:,2) == 3);
immobileOnePre = sum(diffPreMat(:,3) == 0);
confinedOnePre = sum(diffPreMat(:,3) == 1);
brownianOnePre = sum(diffPreMat(:,3) == 2);
directedOnePre = sum(diffPreMat(:,3) == 3);

asymmetricPost = sum(diffPostMat(:,1));
notAsymmetricPost = nInside-asymmetricPost;
immobileFullPost = sum(diffPostMat(:,2) == 0);
confinedFullPost = sum(diffPostMat(:,2) == 1);
brownianFullPost = sum(diffPostMat(:,2) == 2);
directedFullPost = sum(diffPostMat(:,2) == 3);
immobileOnePost = sum(diffPostMat(:,3) == 0);
confinedOnePost = sum(diffPostMat(:,3) == 1);
brownianOnePost = sum(diffPostMat(:,3) == 2);
directedOnePost = sum(diffPostMat(:,3) == 3);

out.nInside = nInside;
out.nOutside = nOutside;
out.totalLengthInside = totalLengthInside;
out.totalLengthOutside = totalLengthOutside;
out.meanLengthInside = meanLengthInside;
out.meanLengthOutside = meanLengthOutside;
out.longestInside = longestInside;
out.longestOutside = longestOutside;
out.nMergeInside = nMergeInside;
out.nMergeOutside = nMergeOutside;
out.nSplitInside =  nSplitInside;
out.nSplitOutside = nSplitOutside;
out.asymmetricInside = asymmetricInside;
out.asymmetricOutside = asymmetricOutside;
out.notAsymmetricInside = notAsymmetricInside;
out.notAsymmetricOutside = notAsymmetricOutside;
out.immobileFullInside = immobileFullInside;
out.immobileFullOutside = immobileFullOutside;
out.confinedFullInside = confinedFullInside;
out.confinedFullOutside = confinedFullOutside;
out.brownianFullInside = brownianFullInside;
out.brownianFullOutside = brownianFullOutside;
out.directedFullInside = directedFullInside;
out.directedFullOutside = directedFullOutside;
out.immobileOneInside = immobileOneInside;
out.immobileOneOutside = immobileOneOutside;
out.confinedOneInside = confinedOneInside;
out.confinedOneOutside = confinedOneOutside;
out.brownianOneInside = brownianOneInside;
out.brownianOneOutside = brownianOneOutside;
out.directedOneInside = directedOneInside;
out.directedOneOutside = directedOneOutside;
out.mssSlopeFullInside = mssSlopeFullInside;
out.mssSlopeFullOutside = mssSlopeFullOutside;
out.genDiffCoefFullInside = genDiffCoefFullInside;
out.genDiffCoeffFullOutside = genDiffCoefFullOutside;
out.scalingPowerFullInside = scalingPowerFullInside;
out.scalingPowerFullOutside = scalingPowerFullOutside;
out.normDiffCoefFullInside = normDiffCoefFullInside;
out.normDiffCoefFullOutside = normDiffCoefFullOutside;
out.mssSlopeOneInside = mssSlopeOneInside;
out.mssSlopeOneOutside = mssSlopeOneOutside;
out.genDiffCoefOneInside = genDiffCoefOneInside;
out.genDiffCoeffOneOutside = genDiffCoefOneOutside;
out.scalingPowerOneInside = scalingPowerOneInside;
out.scalingPowerOneOutside = scalingPowerOneOutside;
out.normDiffCoefOneInside = normDiffCoefOneInside;
out.normDiffCoefOneOutside = normDiffCoefOneOutside;
out.confRadiusInside = confRadiusInside;
out.confRadiusOutside = confRadiusOutside;
out.asymmetricPre = asymmetricPre;
out.asymmetricOutside = asymmetricOutside;
out.notAsymmetricPre = notAsymmetricPre;
out.notAsymmetricOutside = notAsymmetricOutside;
out.immobileFullPre = immobileFullPre;
out.immobileFullPost = immobileFullPost;
out.confinedFullPre = confinedFullPre;
out.confinedFullPost = confinedFullPost;
out.brownianFullPre = brownianFullPre;
out.brownianFullPost = brownianFullPost;
out.directedFullPre = directedFullPre;
out.directedFullPost = directedFullPost;
out.immobileOnePre = immobileOnePre;
out.immobileOnePost = immobileOnePost;
out.confinedOnePre = confinedOnePre;
out.confinedOnePost = confinedOnePost;
out.brownianOnePre = brownianOnePre;
out.brownianOnePost = brownianOnePost;
out.directedOnePre = directedOnePre;
out.directedOnePost = directedOnePost;


end

