function [result, sortResult] = partitionCoef(CML)
%PARTITIONCOEF determines the partition coefficient from partition fraction and control partition fraction from randomized control
%
%INPUT
%   CML     : Combined Movie List or Movie List
%
%OUTPUT
%   result
%       .partCoef   : structure array partition coefficient of MD's. Each
%                     element stores information about a single MD.
%           .immobile
%           .confined
%           .free
%           .directeed
%           .undetermined
%       *All following structures contain above structure elements
%       .coefMean   : structure of weighted mean partition coefficient of MD's
%       .MDWeight   : structure array of weight used to calculate above weighted mean.
%                     The weight is equal to the total number of times all
%                     relevant racks were observed. ie). if we see 2 immobile
%                     tracks: one for 5 frames and the other for 20 frames.
%                     Then, the weight of that MD's imbile component is 25
%                     frames.
%       .meanPart   : structure array of weighted mean of partitioning fraction
%       .meanCont   : structure array of weighted mean of partitioning fraction
%                     of randomized ocntrol.
%       .locFreq    : structure array of localization frequency of each MD.
%                     Similar to k_on. nLocE = k_on * [weight] * [1-meanPart] * [maskboundary]
%                     estimating maskboundary = [meanCont] * [1-meanCont]
%       .delocFreq  : delocalization frequency
%                     Similar to k_off. nDelocE = k_off * [weight] * [meanPart] * [maskboundary]
%       .locFreqMean
%       .delocFreqMean
%       .eqCond     : (equilibrium condition) looks at how close the track
%                     partitioning meets the equilibrium condition. The
%                     equilibrium condition is nLocE - nDelocE == 0. Here I
%                     define eqCond = (nLocE - nDelocE) / weight / [maskboundary].
%
%   sortResult      :result of the analysis
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
%Tae H Kim, August 2015   

%% Initialize
if isa(CML, 'CombinedMovieList')
    if isempty(CML.movieLists_)
        CML.sanityCheck();
    end
    MDs = {};
    nML = numel(CML.movieLists_);
    for iML = 1:nML
        MDs = [MDs CML.movieLists_(iML).movies_]; %#ok<AGROW>
    end
elseif isa(CML, 'MovieList')
    if isempty(CML.movies_)
        pLength = fprintf('Loading MovieList\n');
        evalc('CML.sanityCheck()');
        fprintf(repmat('\b', 1, pLength));
        fprintf('MovieList loaded\n')
    end
    MDs = CML.movies_;
end
nMD = numel(MDs);
FN = {'immobile', 'confined', 'free', 'directed', 'undetermined'};
nFN = numel(FN);

%% Analyze: sorting based on diffusion type
%progress Text
progressTextMultiple('Sorting partition fraction', nMD);
%call analysis
sortResult = cellfun(@callDiffMode, MDs, 'UniformOutput', false);
sortResult = [sortResult{:}];

%% MDWeight
%preallocate
MDWeight(nMD) = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
%get weight
for iMD = 1:nMD
    for iFN = 1:nFN
        MDWeight(iMD).(FN{iFN}) = sum(sortResult(iMD).weightPart.(FN{iFN}));
    end
end

%% meanPart and meanCont
meanPart = [sortResult.wMeanPart];
meanCont = [sortResult.wMeanCont];

%% localization and delocalization events
nLocE = [sortResult.nLocE];
nDelocE = [sortResult.nDelocE];

%% Partition Coefficient
%preallocate
partCoef(nMD) = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
%
for iMD = 1:nMD
    for iFN = 1:nFN
        partCoef(iMD).(FN{iFN}) = getPartPref(meanPart(iMD).(FN{iFN}), meanCont(iMD).(FN{iFN}));
    end
end

for iFN = 1:nFN
    %arithmatic mean
    wSum = 0;
    totW = 0;
    for iMD = 1:nMD
        if ~isnan(partCoef(iMD).(FN{iFN})) && MDWeight(iMD).(FN{iFN}) ~= 0
            wSum = wSum + partCoef(iMD).(FN{iFN}) .* MDWeight(iMD).(FN{iFN});
            totW = totW + MDWeight(iMD).(FN{iFN});
        end
    end
    coefMean.(FN{iFN}) = wSum ./ totW;
    %geometric mean
    totW = 0;
    gSum = 0;
    for iMD = 1:nMD
        if ~isnan(partCoef(iMD).(FN{iFN})) && MDWeight(iMD).(FN{iFN}) ~= 0 && partCoef(iMD).(FN{iFN}) ~= 0
            totW = totW + MDWeight(iMD).(FN{iFN});
            gSum = gSum + log(partCoef(iMD).(FN{iFN})) .* MDWeight(iMD).(FN{iFN});
        end
    end
    coefGMean.(FN{iFN}) = exp(gSum ./ totW);
end

%% Kinetics
%preallocate
locFreq(nMD) = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
delocFreq(nMD) = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
eqCond(nMD) = struct('immobile', [], 'confined', [], 'free', [], 'directed', [], 'undetermined', []);
%frequency determination
for iMD = 1:nMD
    for iFN = 1:nFN
        locFreq(iMD).(FN{iFN}) = nLocE(iMD).(FN{iFN}) ./ MDWeight(iMD).(FN{iFN}) ./ (1 - meanPart(iMD).(FN{iFN})) ./ meanCont(iMD).(FN{iFN}) ./ (1 - meanCont(iMD).(FN{iFN}));
        delocFreq(iMD).(FN{iFN}) = nDelocE(iMD).(FN{iFN}) ./ MDWeight(iMD).(FN{iFN}) ./ meanPart(iMD).(FN{iFN}) ./ meanCont(iMD).(FN{iFN}) ./ (1 - meanCont(iMD).(FN{iFN}));
    end
end
% get mean of locFreq and delocFreq
for iFN = 1:nFN
    locFreqMean.(FN{iFN}) = mean([locFreq.(FN{iFN})]);
    delocFreqMean.(FN{iFN}) = mean([delocFreq.(FN{iFN})]);
end
%eqCond determination
for iMD = 1:nMD
    for iFN = 1:nFN
        eqCond(iMD).(FN{iFN}) = (nLocE(iMD).(FN{iFN}) - nDelocE(iMD).(FN{iFN})) ./ MDWeight(iMD).(FN{iFN}) ./ meanCont(iMD).(FN{iFN}) ./ (1 - meanCont(iMD).(FN{iFN}));
    end
end

%% Combine into single output
result.partCoef = partCoef;
result.coefMean = coefMean;
result.coefGMean = coefGMean;
result.MDWeight = MDWeight;
result.meanPart = meanPart;
result.meanCont = meanCont;
result.locFreq = locFreq;
result.delocFreq = delocFreq;
result.locFreqMean = locFreqMean;
result.delocFreqMean = delocFreqMean;
result.eqCond = eqCond;

%partCoef, coefMean, coefGMean, MDWeight, meanPart, meanCont

end

%% Analyze: sorting based on diffusion type: calls sorting function + progress text
function result = callDiffMode(MD)
result = diffModePartition(MD, false);
progressTextMultiple();
end

%% Structure Concatenate
function s2 = horicatStruct(s)
names = fieldnames(s);
for i = 1:numel(names)
    s2.(names{i}) = [s.(names{i})];
end
end
function s2 = vertcatStruct(s)
names = fieldnames(s);
for i = 1:numel(names)
    s2.(names{i}) = vertcat(s.(names{i}));
end
end

%% Partition Coefficient Calculation
function result = getPartPref(p, c)
result = p ./ c .* (1 - c) ./ (1 - p);
end









