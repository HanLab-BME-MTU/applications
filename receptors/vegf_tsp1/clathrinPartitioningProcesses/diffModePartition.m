function [result] = diffModePartition(MDML, saveData)
%DIFFMODEPARTITION sperates tracks based on diffusion mode and saves patition fraction, and compares it to the control
%
%INPUT
%   ML          : MovieList Object. All MovieData must have
%                 PartitionAnalysisProcess and MotionAnalysisProcess
%   saveData    : (logical) determines whether the data will be used to create
%                 figures and be saved
%
%OUTPUT
%   result      :result of the analysis
%       .partFrac       : (structure) contains partition fraction of tracks
%                     organized based on diffusion mode
%           .immobile       : (array) contains partition fractions of tracks
%                             with immobile diffusion profile
%           .confined       : (array) contains partition fractions of tracks
%                             with confined diffusion profile
%           .free           : (array) contains partition fractions of tracks
%                             with free diffusion profile
%           .directed       : (array) contains partition fractions of tracks
%                             with directional diffusion profile
%           .undetermined   : (array) contains partition fractions of tracks
%                             with indeterminate diffusion profile
%       .meanPart       : average values of partition fraction
%           .immobile       : average
%           .confined       : average
%           .free           : average
%           .directed       : average
%           .undetermined   : average
%   *Following structure elements contain strucutre elements 'immobile',
%   'confined', 'free', 'directed', and 'undetermined' as above.
%       .weightPart     : number of frames the track was observed
%       .wMeanPart      : weighted average
%       .partCont       : control group eq of partFrac
%       .meanCont       : control group eq of meanPart
%       .weightCont     : control group eq of weightPart
%       .wMeanPart      : control group eq of wMeanPart
%       .KS             : p value obtained from two sample Kolmogorov-
%                         Smirnov test between data and randomized control.
%                         Smaller p means the null hypothesis (two data
%                         sets are same) is rejected.
%       .nLocE          : number of localization event (of entire MD or ML)
%       .nDelocE        : number of delocalization event
%
%Tae H Kim, July 2015

%% Initialization
if saveData
    path = uigetdir('Select a directory for saving figures');
    name = input('What is the title?\n', 's');
    %creates directory if doesn't exist
    if ~exist(path, 'dir')
        mkdir(path);
    end
end
%Check input
assert(isa(MDML, 'MovieList') || isa(MDML, 'MovieData'));
%nMD = numel(MDML.movieDataFile_);
%Check if ML is loaded. If not, load.
if isa(MDML, 'MovieList') && isempty(MDML.movies_)
    pLength = fprintf('Loading MovieList\n');
    evalc('MDML.sanityCheck()');
    fprintf(repmat('\b', 1, pLength));
    fprintf('MovieList loaded\n')
end
if nargin < 2
    saveData = false;
end
%struct field
FN = {'immobile', 'confined', 'free', 'directed', 'undetermined'};
nFN = 5;

%% Analysis
%progress display
%progressTextMultiple('analyzing MovieList', nMD);
%deals out MovieData
if isa(MDML, 'MovieList')
    [partFrac_, weightPart_, partCont_, weightCont_] = cellfun(@analyzeMD, MDML.movies_, 'UniformOutput', false);
    nLocE = nan;
    nDelocE = nan;
elseif isa(MDML, 'MovieData')
    [partFrac_{1}, weightPart_{1}, partCont_{1}, weightCont_{1}, nLocE, nDelocE] = analyzeMD(MDML);
end
%combine output (partition fraction: experimental)
partFrac_ = [partFrac_{:}];
partFrac.immobile = horzcat(partFrac_.immobile);
partFrac.confined = horzcat(partFrac_.confined);
partFrac.free = horzcat(partFrac_.free);
partFrac.directed = horzcat(partFrac_.directed);
partFrac.undetermined = horzcat(partFrac_.undetermined);
%combine output (partition weight: experimental)
weightPart_ = [weightPart_{:}];
weightPart.immobile = horzcat(weightPart_.immobile);
weightPart.confined = horzcat(weightPart_.confined);
weightPart.free = horzcat(weightPart_.free);
weightPart.directed = horzcat(weightPart_.directed);
weightPart.undetermined = horzcat(weightPart_.undetermined);
%combine output (partition fraction: randomized control)
partCont_ = [partCont_{:}];
partCont.immobile = horzcat(partCont_.immobile);
partCont.confined = horzcat(partCont_.confined);
partCont.free = horzcat(partCont_.free);
partCont.directed = horzcat(partCont_.directed);
partCont.undetermined = horzcat(partCont_.undetermined);
%combine output (partition weight: randomized control)
weightCont_ = [weightCont_{:}];
weightCont.immobile = horzcat(weightCont_.immobile);
weightCont.confined = horzcat(weightCont_.confined);
weightCont.free = horzcat(weightCont_.free);
weightCont.directed = horzcat(weightCont_.directed);
weightCont.undetermined = horzcat(weightCont_.undetermined);

%% Get unweighted mean
%experimental
meanPart.immobile = mean(partFrac.immobile, 'omitnan');
meanPart.confined = mean(partFrac.confined, 'omitnan');
meanPart.free = mean(partFrac.free, 'omitnan');
meanPart.directed = mean(partFrac.directed, 'omitnan');
meanPart.undetermined = mean(partFrac.undetermined, 'omitnan');
%randomized control
meanCont.immobile = mean(partCont.immobile, 'omitnan');
meanCont.confined = mean(partCont.confined, 'omitnan');
meanCont.free = mean(partCont.free, 'omitnan');
meanCont.directed = mean(partCont.directed, 'omitnan');
meanCont.undetermined = mean(partCont.undetermined, 'omitnan');

%% Get Weighted Mean
%experimental
wMeanPart.immobile = sum(weightPart.immobile .* partFrac.immobile, 'omitnan') ./ sum(weightPart.immobile);
wMeanPart.confined = sum(weightPart.confined .* partFrac.confined, 'omitnan') ./ sum(weightPart.confined);
wMeanPart.free = sum(weightPart.free .* partFrac.free, 'omitnan') ./ sum(weightPart.free);
wMeanPart.directed = sum(weightPart.directed .* partFrac.directed, 'omitnan') ./ sum(weightPart.directed);
wMeanPart.undetermined = sum(weightPart.undetermined .* partFrac.undetermined, 'omitnan') ./ sum(weightPart.undetermined);
%control
wMeanCont.immobile = sum(weightCont.immobile .* partCont.immobile, 'omitnan') ./ sum(weightCont.immobile);
wMeanCont.confined = sum(weightCont.confined .* partCont.confined, 'omitnan') ./ sum(weightCont.confined);
wMeanCont.free = sum(weightCont.free .* partCont.free, 'omitnan') ./ sum(weightCont.free);
wMeanCont.directed = sum(weightCont.directed .* partCont.directed, 'omitnan') ./ sum(weightCont.directed);
wMeanCont.undetermined = sum(weightCont.undetermined .* partCont.undetermined, 'omitnan') ./ sum(weightCont.undetermined);

%% KS-test
try
    [~, KS.immobile] = kstest2(partFrac.immobile, partCont.immobile);
catch
    KS.immobile = nan;
end
try
    [~, KS.confined] = kstest2(partFrac.confined, partCont.confined);
catch
    KS.confined = nan;
end
try
    [~, KS.free] = kstest2(partFrac.free, partCont.free);
catch
    KS.free = nan;
end
try
    [~, KS.directed] = kstest2(partFrac.directed, partCont.directed);
catch
    KS.directed = nan;
end
try
    [~, KS.undetermined] = kstest2(partFrac.undetermined, partCont.undetermined);
catch
    KS.undetermined = nan;
end

%% Combine Data
%Combined all result into one structure
result.partFrac = partFrac;
result.weightPart = weightPart;
result.meanPart = meanPart;
result.wMeanPart = wMeanPart;
result.partCont = partCont;
result.weightCont = weightCont;
result.meanCont = meanCont;
result.wMeanCont = wMeanCont;
result.KS = KS;
result.nLocE = nLocE;
result.nDelocE = nDelocE;

%% Plot and save
if saveData
    %Plot and save figures -------------------------------------------------------------------------------------------------------------------------------
    %immobile
    figObj = figure;
    histogram(partFrac.immobile, 0:.05:1, 'Normalization', 'probability');
    hold on;
    histogram(partCont.immobile, 0:.05:1, 'Normalization', 'probability');
    title(['Partition Fraction Histogram: ' name ' immobile']);
    legend('data', 'randomized control');
    savefig(figObj, [path filesep name ' immobile.fig']);
    close(figObj);
    %confined
    figObj = figure;
    histogram(partFrac.confined, 0:.05:1, 'Normalization', 'probability');
    hold on;
    histogram(partCont.confined, 0:.05:1, 'Normalization', 'probability');
    title(['Partition Fraction Histogram: ' name ' confined']);
    legend('data', 'randomized control');
    savefig(figObj, [path filesep name ' confined.fig']);
    close(figObj);
    %free
    figObj = figure;
    histogram(partFrac.free, 0:.05:1, 'Normalization', 'probability');
    hold on;
    histogram(partCont.free, 0:.05:1, 'Normalization', 'probability');
    title(['Partition Fraction Histogram: ' name ' free']);
    legend('data', 'randomized control');
    savefig(figObj, [path filesep name ' free.fig']);
    close(figObj);
    %directed
    figObj = figure;
    histogram(partFrac.directed, 0:.05:1, 'Normalization', 'probability');
    hold on;
    histogram(partCont.directed, 0:.05:1, 'Normalization', 'probability');
    title(['Partition Fraction Histogram: ' name ' directed']);
    legend('data', 'randomized control');
    savefig(figObj, [path filesep name ' directed.fig']);
    close(figObj);
    %undetermined
    figObj = figure;
    histogram(partFrac.undetermined, 0:.05:1, 'Normalization', 'probability');
    hold on;
    histogram(partCont.undetermined, 0:.05:1, 'Normalization', 'probability');
    title(['Partition Fraction Histogram: ' name ' undetermined']);
    legend('data', 'randomized control');
    savefig(figObj, [path filesep name ' undetermined.fig']);
    close(figObj);

    %save data
    save([path filesep name '_data.mat'], 'result');
end

end
%% Local functions
% analyzes MD
function [partFrac, weightPart, partCont, weightCont, nLocE, nDelocE] = analyzeMD(MD)
%struct field
FN = {'immobile', 'confined', 'free', 'directed', 'undetermined'};
nFN = 5;
%load PartitionAnalysisProcess
processIndx = MD.getProcessIndex('PartitionAnalysisProcess');
load(MD.processes_{processIndx}.outFilePaths_{1});
%load MotionAnalysisProcess
processIndx = MD.getProcessIndex('MotionAnalysisProcess');
load(MD.processes_{processIndx}.outFilePaths_{1});
track = arrayfun(@(x) x.classification(:,2), tracks, 'UniformOutput', false);
%analyze
nTrack = numel(track);
[partFrac, weightPart, nLocE, nDelocE] = sortTracks(partitionResult, nTrack, track);
[partCont, weightCont] = cellfun(@(x) sortTracks(x, nTrack, track), partitionControl, 'UniformOutput', false);
%consolidate partCont and weightCont
partCont = [partCont{:}];
weightCont = [weightCont{:}];
partCont = concatStruct(partCont);
weightCont = concatStruct(weightCont);
%progressText
%progressTextMultiple();
end

%% Concatenates structures
function result = concatStruct(input)
names = fieldnames(input);
nField = numel(names);
for iField = 1:nField
    result.(names{iField}) = [input.(names{iField})];
end
end

%% Analyzes each set of tracks
function [fraction, weight, nLocE, nDelocE] = sortTracks(partition, nTrack, track)
%initialize
fraction = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
weight = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
nLocE = struct('immobile', 0, 'confined', 0, 'free', 0, 'directed', 0, 'undetermined', 0);
nDelocE = struct('immobile', 0, 'confined', 0, 'free', 0, 'directed', 0, 'undetermined', 0);
for iTrack = 1:nTrack
    nSubTrack = numel(partition(iTrack).partitionFrac);
    for iSubTrack = 1:nSubTrack
        diffMode = track{iTrack}(iSubTrack);
        if isnan(diffMode)
            fraction.undetermined(end+1) = partition(iTrack).partitionFrac(iSubTrack);
            weight.undetermined(end+1) = partition(iTrack).nFramesTot(iSubTrack);
            nLocE.undetermined = nLocE.undetermined + partition(iTrack).nLocEvent(iSubTrack);
            nDelocE.undetermined = nDelocE.undetermined + partition(iTrack).nDelocEvent(iSubTrack);
        elseif diffMode == 0
            fraction.immobile(end+1) = partition(iTrack).partitionFrac(iSubTrack);
            weight.immobile(end+1) = partition(iTrack).nFramesTot(iSubTrack);
            nLocE.immobile = nLocE.immobile + partition(iTrack).nLocEvent(iSubTrack);
            nDelocE.immobile = nDelocE.immobile + partition(iTrack).nDelocEvent(iSubTrack);
        elseif diffMode == 1
            fraction.confined(end+1) = partition(iTrack).partitionFrac(iSubTrack);
            weight.confined(end+1) = partition(iTrack).nFramesTot(iSubTrack);
            nLocE.confined = nLocE.confined + partition(iTrack).nLocEvent(iSubTrack);
            nDelocE.confined = nDelocE.confined + partition(iTrack).nDelocEvent(iSubTrack);
        elseif diffMode == 2
            fraction.free(end+1) = partition(iTrack).partitionFrac(iSubTrack);
            weight.free(end+1) = partition(iTrack).nFramesTot(iSubTrack);
            nLocE.free = nLocE.free + partition(iTrack).nLocEvent(iSubTrack);
            nDelocE.free = nDelocE.free + partition(iTrack).nDelocEvent(iSubTrack);
        elseif diffMode == 3
            fraction.directed(end+1) = partition(iTrack).partitionFrac(iSubTrack);
            weight.directed(end+1) = partition(iTrack).nFramesTot(iSubTrack);
            nLocE.directed = nLocE.directed + partition(iTrack).nLocEvent(iSubTrack);
            nDelocE.directed = nDelocE.directed + partition(iTrack).nDelocEvent(iSubTrack);
        end
    end
end
end












