function [partFrac, meanPart] = diffModePartition(ML, saveData)
%DIFFMODEPARTITION sperates tracks based on diffusion mode and saves patition fraction
%
%INPUT
%   ML          : MovieList Object. All MovieData must have
%                 PartitionAnalysisProcess and MotionAnalysisProcess
%   saveData    : (logical) determines whether the data will be used to create
%                 figures and be saved
%
%OUTPUT
%   partFrac    : (structure) contains partition fraction of tracks
%                 organized based on diffusion mode
%       .immobile       : (array) contains partition fractions of tracks
%                         with immobile diffusion profile
%       .confined       : (array) contains partition fractions of tracks
%                         with confined diffusion profile
%       .free           : (array) contains partition fractions of tracks
%                         with free diffusion profile
%       .directed       : (array) contains partition fractions of tracks
%                         with directional diffusion profile
%       .undetermined   : (array) contains partition fractions of tracks
%                         with indeterminate diffusion profile
%   meanPart    : average values of partition fraction
%       .immobile       : average
%       .confined       : average
%       .free           : average
%       .directed       : average
%       .undetermined   : average
%
%Tae H Kim, July 2015

%% Initialization
if saveData
    path = uigetdir('Select a directory for saving figures');
    name = input('What is the title?\n', 's');
end
%Check input
assert(isa(ML, 'MovieList'));
nMD = numel(ML.movieDataFile_);
%Check if ML is loaded. If not, load.
if isempty(ML.movies_)
    pLength = fprintf('Loading MovieList\n');
    evalc('ML.sanityCheck()');
    fprintf(repmat('\b', 1, pLength));
    fprintf('MovieList loaded\n')
end
if nargin < 2
    saveData = false;
end

%% Analysis
%progress display
%multipleProgressText('analyzing MovieList', nMD);
%deals out MovieData
partFrac_ = cellfun(@analyzeMD, ML.movies_, 'UniformOutput', false);
%combine output
partFrac_ = [partFrac_{:}];
partFrac.immobile = horzcat(partFrac_.immobile);
partFrac.confined = horzcat(partFrac_.confined);
partFrac.free = horzcat(partFrac_.free);
partFrac.directed = horzcat(partFrac_.directed);
partFrac.undetermined = horzcat(partFrac_.undetermined);

%% Further analysis and save
meanPart.immobile = mean(partFrac.immobile, 'omitnan');
meanPart.confined = mean(partFrac.confined, 'omitnan');
meanPart.free = mean(partFrac.free, 'omitnan');
meanPart.directed = mean(partFrac.directed, 'omitnan');
meanPart.undetermined = mean(partFrac.undetermined, 'omitnan');
if saveData
    %immobile
    figObj = figure;
    histogram(partFrac.immobile, 0:.05:1);
    title(['Partition Fraction Histogram: ' name ' immobile']);
    savefig(figObj, [path filesep name ' immobile.fig']);
    close(figObj);
    %confined
    figObj = figure;
    histogram(partFrac.confined, 0:.05:1);
    title(['Partition Fraction Histogram: ' name ' confined']);
    savefig(figObj, [path filesep name ' confined.fig']);
    close(figObj);
    %free
    figObj = figure;
    histogram(partFrac.free, 0:.05:1);
    title(['Partition Fraction Histogram: ' name ' free']);
    savefig(figObj, [path filesep name ' free.fig']);
    close(figObj);
    %directed
    figObj = figure;
    histogram(partFrac.directed, 0:.05:1);
    title(['Partition Fraction Histogram: ' name ' directed']);
    savefig(figObj, [path filesep name ' directed.fig']);
    close(figObj);
    %undetermined
    figObj = figure;
    histogram(partFrac.undetermined, 0:.05:1);
    title(['Partition Fraction Histogram: ' name ' undetermined']);
    savefig(figObj, [path filesep name ' undetermined.fig']);
    close(figObj);
    %save data
    save([path filesep name '_data.mat'], 'partFrac', 'meanPart');
end

end
%% Local functions
% analyzes MD
function [partFrac] = analyzeMD(MD)
%initialize
partFrac = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
%load PartitionAnalysisProcess
processIndx = MD.getProcessIndex('PartitionAnalysisProcess');
load(MD.processes_{processIndx}.outFilePaths_{1});
%load MotionAnalysisProcess
processIndx = MD.getProcessIndex('MotionAnalysisProcess');
load(MD.processes_{processIndx}.outFilePaths_{1});
nTrack = numel(tracks);
for iTrack = 1:nTrack
    diffMode = tracks(iTrack).classification(1, 2);
    if isnan(diffMode)
        partFrac.undetermined(end+1) = partitionResult(iTrack).partitionFrac;
    elseif diffMode == 0
        partFrac.immobile(end+1) = partitionResult(iTrack).partitionFrac;
    elseif diffMode == 1
        partFrac.confined(end+1) = partitionResult(iTrack).partitionFrac;
    elseif diffMode == 2
        partFrac.free(end+1) = partitionResult(iTrack).partitionFrac;
    elseif diffMode == 3
        partFrac.directed(end+1) = partitionResult(iTrack).partitionFrac;
    end
end
%progressText
%multipleProgressText();
end












