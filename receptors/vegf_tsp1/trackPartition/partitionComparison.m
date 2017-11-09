function out = partitionComparison(varargin)
%PARTITIONCOMPARISON Aggregate diffusion analysis and merge/split
%information from all partitioned tracks in a ML or CML and plot the
%results. If run with no arguments, the function will prompt the user to
%select ML/CMLs, input condition names, and select an output directory for
%saving plots and aggregated data.
%   
% Inputs: string name of condition 1, path to CML1, string name of condition 2,
% path to CML2, ...
% User will be prompted to select an output directory.

% Collect input MovieLists into one cell array
% Each cell of this array contains a cell array of information from the
% movies under each condition
CML = {};
conditions = {};
if (nargin == 0) || (mod(nargin,2) ~= 0)
    [file,path] = uigetfile('*.mat','Select a CML for condition 1');
    if file == 0
        error('No CML selected')
    end
    cond = inputdlg('Name condition 1:');
    CML = [CML,{[path,file]}];
    conditions = [conditions,cond];
    while file ~= 0
        condNum = numel(CML);
        [file,path] = uigetfile('*.mat',sprintf('Select a CML for condition %g',condNum+1));
        if file ~= 0
            cond = inputdlg(sprintf('Name condition %g:',condNum+1));
            CML = [CML,{[path,file]}];
            conditions = [conditions,cond];
        end
    end     
else
    CML = varargin(2*(1:nargin/2));
    conditions = varargin(2*(1:nargin/2)-1);
end

outDir = uigetdir('','Select directory for saving plots');

MDInfo = cell(1,numel(CML));
for i = 1:numel(CML)
        file = load(CML{i});
        if isfield(file,'CML')
            MDInfo{i} = loadCML(file.CML);
        elseif isfield(file,'ML')
            MDInfo{i} = loadML(file.ML);
        else
            error('Input %g is not a valid ML or CML path',i)
        end
end

out = cell(1,numel(CML));
for i = 1:numel(CML)
    [out{i},insidePrePost{i}] = parcellfun_progress(@getInfo,MDInfo{i}','UniformOutput',false);
    out{i} = cell2mat(out{i});
end

save([outDir,filesep,'partitionComparison.mat'],'out','insidePrePost')
plotComparisonInfo(out,conditions,outDir);
end

function MDInfo = loadML(ML)
% Load track partition results from all MDs in a ML into a cell array

nMD = numel(ML.movieDataFile_);
MDInfo = cell(nMD,1);
valid = ones(nMD,1);
for i = 1:nMD
    file = load(ML.movieDataFile_{i});
    if file.MD.getProcessIndex('TrackPartitionProcess') > 0
        process = file.MD.processes_{file.MD.getProcessIndex('TrackPartitionProcess')};
        processOut = load(process.outFilePaths_{1});
        MDInfo{i}.tracks = processOut.tracksPart;
        MDInfo{i}.diffAnalysisRes = processOut.diffAnalysisRes;
        MDInfo{i}.info = processOut.info;
    else
        % If MD does not have TrackPartitionProcess, skip it
        valid(i) = 0;
    end
end
% Use only MDs with TrackPartitionProcess
MDInfo = MDInfo(find(valid));
% Make sure it's a column
MDInfo = MDInfo(:);
end

function MDInfo = loadCML(CML)
% Load track partition results from all MDs in a CML into a cell array
nML = numel(CML.movieListDirectory_);
MDInfo = {};
for i = 1:nML
    file = load(CML.movieListDirectory_{i});
    trackInfoTemp = loadML(file.ML);
    MDInfo = [MDInfo;trackInfoTemp];
end
end

function [out,insidePrePost] = getInfo(MDInfo)
% Aggregate the diffusion analysis info from all the MDs

isInside = arrayfun(@(x) x.isInside,MDInfo.tracks);
isOutside = ~isInside;
tracksInside = MDInfo.tracks(isInside);
tracksOutside = MDInfo.tracks(isOutside);
diffInside = MDInfo.diffAnalysisRes(isInside);
diffOutside = MDInfo.diffAnalysisRes(isOutside);

seqInside = arrayfun(@(x) x.seqOfEvents,tracksInside,'UniformOutput',false);
seqInsideMat = cell2mat(seqInside(:));
seqOutside = arrayfun(@(x) x.seqOfEvents,tracksOutside,'UniformOutput',false);
seqOutsideMat = cell2mat(seqOutside(:));

% Count all tracks (including those inside compound tracks)
nInside = size(seqInsideMat,1)/2;
nOutside = size(seqOutsideMat,1)/2;
nTracks = nInside+nOutside;

% In case there are no inside tracks, create an empty seqOfEvents
if isempty(seqInsideMat)
    seqInsideMat = [0 0 0 0];
    nInside = 0;
end

% Count merges and splits
mergeInside = sum((seqInsideMat(:,2) == 2)&(~isnan(seqInsideMat(:,4)))); 
splitInside = sum((seqInsideMat(:,2) == 1)&(~isnan(seqInsideMat(:,4))));
mergeOutside = sum((seqOutsideMat(:,2) == 2)&(~isnan(seqOutsideMat(:,4)))); 
splitOutside = sum((seqOutsideMat(:,2) == 1)&(~isnan(seqOutsideMat(:,4))));
fracMergeInside = mergeInside/nInside;
fracSplitInside = splitInside/nInside;
fracMergeOutside = mergeOutside/nOutside;
fracSplitOutside = splitOutside/nOutside;

% Find mean length of tracks
lengthInside = zeros(nInside,1);
for i = 1:numel(tracksInside)
    lengthInside(i) = numel(tracksInside(i).tracksFeatIndxCG);
end
totalLengthInside = sum(lengthInside);
meanLengthInside = mean(lengthInside);
longestInside = max(lengthInside);

lengthOutside = zeros(nOutside,1);
for i = 1:numel(tracksOutside)
    lengthOutside(i) = numel(tracksOutside(i).tracksFeatIndxCG);
end
totalLengthOutside = sum(lengthOutside);
meanLengthOutside = mean(lengthOutside);
longestOutside = max(lengthOutside);

% Concatenate the diffusion analysis results into big matrices
diffInsideCell = struct2cell(diffInside)';
diffInsideMat = cell2mat(diffInsideCell(:,1));
diffOutsideCell = struct2cell(diffOutside)';
diffOutsideMat = cell2mat(diffOutsideCell(:,1));

if numel(diffInsideMat) == 0
    diffInsideMat = [0,0,0];
end
if numel(diffOutsideMat) == 0
    diffOutsideMat = [0,0,0];
end

% Find the proportion of tracks in each classification, omitting tracks
% with NaN values (unclassified)
asymmetricInside = mean(nonan(diffInsideMat(:,1)) == 1);
notAsymmetricInside = mean(nonan(diffInsideMat(:,1)) == 0);
immobileFullInside = mean(nonan(diffInsideMat(:,2)) == 0);
confinedFullInside = mean(nonan(diffInsideMat(:,2)) == 1);
brownianFullInside = mean(nonan(diffInsideMat(:,2)) == 2);
directedFullInside = mean(nonan(diffInsideMat(:,2)) == 3);
immobileOneInside = mean(nonan(diffInsideMat(:,3)) == 0);
confinedOneInside = mean(nonan(diffInsideMat(:,3)) == 1);
brownianOneInside = mean(nonan(diffInsideMat(:,3)) == 2);
directedOneInside = mean(nonan(diffInsideMat(:,3)) == 3);
% Number of tracks that were actually able to be classified (don't have NaN
% values)
asymmetricInsideN = numel(nonan(diffInsideMat(:,1)) == 1);
notAsymmetricInsideN = numel(nonan(diffInsideMat(:,1)));
immobileFullInsideN = numel(nonan(diffInsideMat(:,2)));
confinedFullInsideN = numel(nonan(diffInsideMat(:,2)));
brownianFullInsideN = numel(nonan(diffInsideMat(:,2)));
directedFullInsideN = numel(nonan(diffInsideMat(:,2)));
immobileOneInsideN = numel(nonan(diffInsideMat(:,3)));
confinedOneInsideN = numel(nonan(diffInsideMat(:,3)));
brownianOneInsideN = numel(nonan(diffInsideMat(:,3)));
directedOneInsideN = numel(nonan(diffInsideMat(:,3)));

asymmetricOutside = mean(nonan(diffOutsideMat(:,1)) == 1);
notAsymmetricOutside = mean(nonan(diffOutsideMat(:,1)) == 0);
immobileFullOutside = mean(nonan(diffOutsideMat(:,2)) == 0);
confinedFullOutside = mean(nonan(diffOutsideMat(:,2)) == 1);
brownianFullOutside = mean(nonan(diffOutsideMat(:,2)) == 2);
directedFullOutside = mean(nonan(diffOutsideMat(:,2)) == 3);
immobileOneOutside = mean(nonan(diffOutsideMat(:,3)) == 0);
confinedOneOutside = mean(nonan(diffOutsideMat(:,3)) == 1);
brownianOneOutside = mean(nonan(diffOutsideMat(:,3)) == 2);
directedOneOutside = mean(nonan(diffOutsideMat(:,3)) == 3);
asymmetricOutsideN = numel(nonan(diffOutsideMat(:,1)) == 1);

notAsymmetricOutsideN = numel(nonan(diffOutsideMat(:,1)));
immobileFullOutsideN = numel(nonan(diffOutsideMat(:,2)));
confinedFullOutsideN = numel(nonan(diffOutsideMat(:,2)));
brownianFullOutsideN = numel(nonan(diffOutsideMat(:,2)));
directedFullOutsideN = numel(nonan(diffOutsideMat(:,2)));
immobileOneOutsideN = numel(nonan(diffOutsideMat(:,3)));
confinedOneOutsideN = numel(nonan(diffOutsideMat(:,3)));
brownianOneOutsideN = numel(nonan(diffOutsideMat(:,3)));
directedOneOutsideN = numel(nonan(diffOutsideMat(:,3)));

% Compute means of these MSS results, omitting NaN values
mssSlopeFullInside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.mssSlope,diffInside,'UniformOutput',false)));
genDiffCoefFullInside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.genDiffCoef,diffInside,'UniformOutput',false)));
scalingPowerFullInside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.scalingPower,diffInside,'UniformOutput',false)));
normDiffCoefFullInside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.normDiffCoef,diffInside,'UniformOutput',false)));
mssSlopeOneInside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.mssSlope,diffInside,'UniformOutput',false)));
genDiffCoefOneInside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.genDiffCoef,diffInside,'UniformOutput',false)));
scalingPowerOneInside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.scalingPower,diffInside,'UniformOutput',false)));
normDiffCoefOneInside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.normDiffCoef,diffInside,'UniformOutput',false)));
confRadiusInside = nanmean(cell2mat(arrayfun(@(x) x.confRadInfo.confRadius,diffInside,'UniformOutput',false)));

mssSlopeFullOutside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.mssSlope,diffOutside,'UniformOutput',false)));
genDiffCoefFullOutside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.genDiffCoef,diffOutside,'UniformOutput',false)));
scalingPowerFullOutside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.scalingPower,diffOutside,'UniformOutput',false)));
normDiffCoefFullOutside = nanmean(cell2mat(arrayfun(@(x) x.fullDim.normDiffCoef,diffOutside,'UniformOutput',false)));
mssSlopeOneOutside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.mssSlope,diffOutside,'UniformOutput',false)));
genDiffCoefOneOutside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.genDiffCoef,diffOutside,'UniformOutput',false)));
scalingPowerOneOutside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.scalingPower,diffOutside,'UniformOutput',false)));
normDiffCoefOneOutside = nanmean(cell2mat(arrayfun(@(x) x.oneDim.normDiffCoef,diffOutside,'UniformOutput',false)));
confRadiusOutside = nanmean(cell2mat(arrayfun(@(x) x.confRadInfo.confRadius,diffOutside,'UniformOutput',false)));

% Get tracks that occur pre or post-inside tracks
% First, get origin track numbers (the index of the original,
% un-partitioned tracks within their track struct)
trackOrigins = arrayfun(@(x) x.originCompoundTrack,MDInfo.tracks);
insideInd = find(isInside);
pre = [];
post = [];
% Save as output variable in this struct:
insidePrePost = struct('insideTrack',[],'insideTrackDiff',[],'preTrack',[],...
    'preTrackDiff',[],'postTrack',[],'postTrackDiff',[]);
for i = 1:numel(insideInd)
    ind = insideInd(i); % index of this track within the entire track struct
    origin = MDInfo.tracks(ind).originCompoundTrack;
    siblings = find(trackOrigins == origin); % find indices of other tracks with same origin
    siblingsPre = siblings(siblings < ind); % find sibling tracks before this track
    siblingsPost = siblings(siblings > ind); % find sibling tracks after this track
    
    pre = [pre,siblingsPre(:)'];
    post = [post,siblingsPost(:)'];
    
    insidePrePost(i).insideTrack = MDInfo.tracks(ind);
    insidePrePost(i).insideTrackDiff = MDInfo.diffAnalysisRes(ind);
    insidePrePost(i).preTrack = MDInfo.tracks(siblingsPre);
    insidePrePost(i).preTrackDiff = MDInfo.diffAnalysisRes(siblingsPre);
    insidePrePost(i).postTrack = MDInfo.tracks(siblingsPost);
    insidePrePost(i).postTrackDiff = MDInfo.diffAnalysisRes(siblingsPost);
    
end

% Make sure this is a column
insidePrePost = insidePrePost(:);

diffPre = MDInfo.diffAnalysisRes(pre);
diffPost = MDInfo.diffAnalysisRes(post);

if isempty(diffPre)
    diffPreMat = [NaN,NaN,NaN,NaN];
else
    diffPreCell = struct2cell(diffPre)';
    diffPreMat = cell2mat(diffPreCell(:,1));
end
if isempty(diffPost)
    diffPostMat = [NaN,NaN,NaN,NaN];
else 
    diffPostCell = struct2cell(diffPost)';
    diffPostMat = cell2mat(diffPostCell(:,1));
end
if numel(diffPreMat) == 0
    diffPreMat = [0,0,0];
end
if numel(diffPostMat) == 0
    diffPostMat = [0,0,0];
end

asymmetricPre = mean(nonan(diffPreMat(:,1)) == 1);
notAsymmetricPre = mean(nonan(diffPreMat(:,1)) == 0);
immobileFullPre = mean(nonan(diffPreMat(:,2)) == 0);
confinedFullPre = mean(nonan(diffPreMat(:,2)) == 1);
brownianFullPre = mean(nonan(diffPreMat(:,2)) == 2);
directedFullPre = mean(nonan(diffPreMat(:,2)) == 3);
immobileOnePre = mean(nonan(diffPreMat(:,3)) == 0);
confinedOnePre = mean(nonan(diffPreMat(:,3)) == 1);
brownianOnePre = mean(nonan(diffPreMat(:,3)) == 2);
directedOnePre = mean(nonan(diffPreMat(:,3)) == 3);
asymmetricPreN = numel(nonan(diffPreMat(:,1)) == 1);
notAsymmetricPreN = numel(nonan(diffPreMat(:,1)));
immobileFullPreN = numel(nonan(diffPreMat(:,2)));
confinedFullPreN = numel(nonan(diffPreMat(:,2)));
brownianFullPreN = numel(nonan(diffPreMat(:,2)));
directedFullPreN = numel(nonan(diffPreMat(:,2)));
immobileOnePreN = numel(nonan(diffPreMat(:,3)));
confinedOnePreN = numel(nonan(diffPreMat(:,3)));
brownianOnePreN = numel(nonan(diffPreMat(:,3)));
directedOnePreN = numel(nonan(diffPreMat(:,3)));

asymmetricPost = mean(nonan(diffPostMat(:,1)) == 1);
notAsymmetricPost = mean(nonan(diffPostMat(:,1)) == 0);
immobileFullPost = mean(nonan(diffPostMat(:,2)) == 0);
confinedFullPost = mean(nonan(diffPostMat(:,2)) == 1);
brownianFullPost = mean(nonan(diffPostMat(:,2)) == 2);
directedFullPost = mean(nonan(diffPostMat(:,2)) == 3);
immobileOnePost = mean(nonan(diffPostMat(:,3)) == 0);
confinedOnePost = mean(nonan(diffPostMat(:,3)) == 1);
brownianOnePost = mean(nonan(diffPostMat(:,3)) == 2);
directedOnePost = mean(nonan(diffPostMat(:,3)) == 3);
asymmetricPostN = numel(nonan(diffPostMat(:,1)) == 1);
notAsymmetricPostN = numel(nonan(diffPostMat(:,1)));
immobileFullPostN = numel(nonan(diffPostMat(:,2)));
confinedFullPostN = numel(nonan(diffPostMat(:,2)));
brownianFullPostN = numel(nonan(diffPostMat(:,2)));
directedFullPostN = numel(nonan(diffPostMat(:,2)));
immobileOnePostN = numel(nonan(diffPostMat(:,3)));
confinedOnePostN = numel(nonan(diffPostMat(:,3)));
brownianOnePostN = numel(nonan(diffPostMat(:,3)));
directedOnePostN = numel(nonan(diffPostMat(:,3)));

out.fracInside = nInside/nTracks;
out.fracOutside = nOutside/nTracks;
out.totalLengthInside = totalLengthInside;
out.totalLengthOutside = totalLengthOutside;
out.meanLengthInside = meanLengthInside;
out.meanLengthOutside = meanLengthOutside;
out.longestInside = longestInside;
out.longestOutside = longestOutside;
out.fracMergeInside = fracMergeInside;
out.fracMergeOutside = fracMergeOutside;
out.fracSplitInside =  fracSplitInside;
out.fracSplitOutside = fracSplitOutside;
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
out.asymmetricPost = asymmetricPost;
out.notAsymmetricPre = notAsymmetricPre;
out.notAsymmetricPost = notAsymmetricPost;
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

out.asymmetricInsideN = asymmetricInsideN;
out.asymmetricOutsideN = asymmetricOutsideN;
out.notAsymmetricInsideN = notAsymmetricInsideN;
out.notAsymmetricOutsideN = notAsymmetricOutsideN;
out.immobileFullInsideN = immobileFullInsideN;
out.immobileFullOutsideN = immobileFullOutsideN;
out.confinedFullInsideN = confinedFullInsideN;
out.confinedFullOutsideN = confinedFullOutsideN;
out.brownianFullInsideN = brownianFullInsideN;
out.brownianFullOutsideN = brownianFullOutsideN;
out.directedFullInsideN = directedFullInsideN;
out.directedFullOutsideN = directedFullOutsideN;
out.immobileOneInsideN = immobileOneInsideN;
out.immobileOneOutsideN = immobileOneOutsideN;
out.confinedOneInsideN = confinedOneInsideN;
out.confinedOneOutsideN = confinedOneOutsideN;
out.brownianOneInsideN = brownianOneInsideN;
out.brownianOneOutsideN = brownianOneOutsideN;
out.directedOneInsideN = directedOneInsideN;
out.directedOneOutsideN = directedOneOutsideN;
out.asymmetricPreN = asymmetricPreN;
out.asymmetricPostN = asymmetricPostN;
out.notAsymmetricPreN = notAsymmetricPreN;
out.notAsymmetricPostN = notAsymmetricPostN;
out.immobileFullPreN = immobileFullPreN;
out.immobileFullPostN = immobileFullPostN;
out.confinedFullPreN = confinedFullPreN;
out.confinedFullPostN = confinedFullPostN;
out.brownianFullPreN = brownianFullPreN;
out.brownianFullPostN = brownianFullPostN;
out.directedFullPreN = directedFullPreN;
out.directedFullPostN = directedFullPostN;
out.immobileOnePreN = immobileOnePreN;
out.immobileOnePostN = immobileOnePostN;
out.confinedOnePreN = confinedOnePreN;
out.confinedOnePostN = confinedOnePostN;
out.brownianOnePreN = brownianOnePreN;
out.brownianOnePostN = brownianOnePostN;
out.directedOnePreN = directedOnePreN;
out.directedOnePostN = directedOnePostN;
end

function out = nonan(vector)
out = vector(~isnan(vector));
end
