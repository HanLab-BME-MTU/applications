function aggregParamMS = receptorAggregParamSimple(tracksAggreg,diffAnalysisRes,...
    timeStep,minTrackLen,probDim,removePotArtifacts,maxClustSize)
%RECEPTORAGGREGPARAM calculates various parameters characterizing receptor aggregation kinetics
%
%SYNOPSIS aggregParamMS = receptorAggregParamSimple(tracksAggreg,diffAnalysisRes,...
%    timeStep,minTrackLen,probDim,removePotArtifacts,maxClustSize)
%
%INPUT  tracksAggreg   : Cell array of output of aggregStateFromCompTracks.
%                        1 entry = 1 movie/simulation.
%       diffAnalysisRes: Cell array of track diffusion analysis results,
%                        i.e. output of trackDiffusionAnalysis1.
%       timeStep       : Time step between frames (s).
%       minTrackLen    : Minimum length of a track to be used in
%                        parameter calculation (frames).
%                        Optional. Default: 5 frames.
%       probDim        : Dimensionality - 2 for 2D, 3 for 3D.
%                        Optional. Default: 2.
%       removePotArtifacts: 1 to remove potentially artifactual merges and
%                        splits, resulting for instance from detection
%                        artifacts, 0 otherwise. 
%                        Optional. Default: 1.
%       maxClustSize   : Maximum cluster size to use for cluster
%                        size distribution analysis. Any clusters
%                        larger than maxClustSize will be lumped into the
%                        bin for maxClustSize.
%                        Optional. Default: maximum cluster size in data.
%
%
%OUTPUT aggregParamMS  : 2D array of mean & std of aggregation parameters.
%                        The rows correspond to:
%                        (1) Rate of merging.
%                        (2) Rate of splitting.
%                        (3) Average merge-to-split time.
%                        (4) Average split-to-merge time.
%                        (5...5+maxClustSize-1) Fraction of clusters of size
%                            1, 2, ..., maxClustSize receptors.
%
%Khuloud Jaqaman, March 2013

%% Input

if nargin < 3
    error('receptorAggregParamSimple: Wrong number of input arguments')
end

%get number of movies
numMovie = length(tracksAggreg);

%assign some default parameters
if nargin < 4 || isempty(minTrackLen)
    minTrackLen = 5;
end
if nargin < 5 || isempty(probDim)
    probDim = 2;
end
if nargin < 6 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

%get maximum cluster size if not input
if nargin < 7 || isempty(maxClustSize)
    aggregMat = cell(numMovie,1);
    for iMovie = 1 : numMovie
        [~,~,~,~,aggregMat{iMovie}] = convStruct2MatIgnoreMS(tracksAggreg{iMovie}.defaultFormatTracks);
    end
    aggregMatAll = vertcat(aggregMat{:});
    maxClustSize = max(aggregMatAll(:));
end

%% Rate of merging and splitting per particle

%get general merge and split statistics
statsGeneral = zeros(numMovie,5);
for iMovie = 1 : numMovie
    statsGeneral(iMovie,:) = calcStatsMS(tracksAggreg{iMovie}.defaultFormatTracks,...
        minTrackLen,probDim,diffAnalysisRes{iMovie},removePotArtifacts);
end

%divide probabilities by time step to get rates
rateMergeAll = statsGeneral(:,4) / timeStep;
rateSplitAll = statsGeneral(:,5) / timeStep;

%calculate mean and std
rateMergeMS = [mean(rateMergeAll) std(rateMergeAll)];
rateSplitMS = [mean(rateSplitAll) std(rateSplitAll)];

%% Distribution of cluster sizes

%calculate normalized distribution of cluster sizes
%this is already mean and std of per movie values
fracClustMS = calcFracClust(tracksAggreg,maxClustSize);

%% Merge-to-split time & split-to-merge time

%calculate merge-to-split time and split-to-merge time
m2sTimeMean = zeros(numMovie,1);
s2mTimeMean = zeros(numMovie,1);
for iMovie = 1 : numMovie
    
    msTimeInfo = calcMergeSplitTimes(tracksAggreg{iMovie}.defaultFormatTracks,...
        minTrackLen,probDim,diffAnalysisRes{iMovie},removePotArtifacts);
    
    %merge-to-split time
    m2sTime = [msTimeInfo.conf.timeMerge2Split; msTimeInfo.brown.timeMerge2Split; ...
        msTimeInfo.lin.timeMerge2Split];
    m2sTimeMean(iMovie) = mean(m2sTime);
    
    %split-to-merge time
    s2mTime = [msTimeInfo.conf.timeSplit2MergeSelf; msTimeInfo.conf.timeSplit2MergeOther; ...
        msTimeInfo.brown.timeSplit2MergeSelf; msTimeInfo.brown.timeSplit2MergeOther; ...
        msTimeInfo.lin.timeSplit2MergeSelf; msTimeInfo.lin.timeSplit2MergeOther];
    s2mTimeMean(iMovie) = mean(s2mTime);
    
end

%calculate mean and std
m2sTimeMS = [mean(m2sTimeMean) std(m2sTimeMean)]*timeStep;
s2mTimeMS = [mean(s2mTimeMean) std(s2mTimeMean)]*timeStep;

%% Output parameter vector

aggregParamMS = [rateMergeMS; rateSplitMS; m2sTimeMS; s2mTimeMS; fracClustMS];

%% ~~~ the end ~~~
