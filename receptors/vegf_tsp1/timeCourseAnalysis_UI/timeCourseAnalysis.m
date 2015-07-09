function [] = timeCourseAnalysis(CMLs, outputDir, varargin)
%TimeCourseAnalysis of CombinedMovieList objects.
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%
%SYNOPSIS function [] = timeCourseAnalysis(CMLs, outputDir, varargin)
%
%INPUT
%   CMLs        : (array of CombinedMovieList) objects or (cell of string)
%                 containingfull paths to CMLs
%   ouputDir    : (string) Directory where figures and analysis result will
%                 be stored
%   varargin    : name_value pairs for analysis parameter
%       'smoothingPara'         : parameter used for smoothing spline fit
%       'channel'               : Channel of MD to be analyzed. Default = 1
%       'doPartitionAnalysis'   : (logical) to analyzes trackPartitioning
%                                 process or not
%
%Tae H Kim, July 2015

%% Input
%convert CMLs if it's cell of strings
if iscellstr(CMLs)
    directory_CML = CMLs;
    CMLs = cellfun(@CombinedMovieList.load, directory_CML, 'UniformOutput', false);
    CMLs = [CMLs{:}];
elseif isa(CMLs, 'CombinedMovieList')
    directory_CML = arrayfun(@(x) x.getFullPath(), CMLs, 'UniformOutput', false); %#ok<NASGU>
end
%input parser
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('CMLs', @(x) isa(x, 'CombinedMovieList'));
ip.addRequired('outputDir', @ischar);
ip.addParameter('smoothingPara', .01, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParameter('channel', 1, @isnumeric);
ip.addParameter('doPartitionAnalysis', false, @(x) isnumeric(x) || islogical(x));
ip.parse(CMLs, outputDir, varargin{:});
smoothingPara = ip.Results.smoothingPara;
analysisPara.smoothingPara = smoothingPara;
channel = round(ip.Results.channel);
%used for Extra Analysis
analysisPara.doPartition = ip.Results.doPartitionAnalysis;
%% Main Time Course Analysis
%For analysis progress display
nCML = numel(CMLs);
%Using resultsIndTimeCourseMod.m to do basic analysis
%and extract time data and align
[summary, time, extra] = arrayfun(@CMLAnalyze, CMLs, 1:nCML, 'UniformOutput', false);
%nested function for above cellfun
%deals with individual CML
    function [CMLSummary, CMLTime, CMLExtra] = CMLAnalyze(CML, iCML)
        alignEvent = CML.analysisPara_.alignEvent;
        %[CMLSummary, CMLTime, CMLExtra] = arrayfun(@(x) MLAnalyze(x, alignEvent), CML.movieLists_, 'UniformOutput', false);
        nML = numel(CML.movieLists_);
        for iML = nML:-1:1
            printLength = fprintf ('Analyzing Movie List %g/%g of CML %g/%g\n', nML-iML+1, nML, iCML, nCML);
            [CMLSummary(iML), CMLTime(iML), CMLExtra(iML)] = MLAnalyze(CML.movieLists_(iML), alignEvent);
            fprintf(repmat('\b',1,printLength));
        end
        CMLSummary = vertcat(CMLSummary{:});
        CMLTime = vertcat(CMLTime{:});
        CMLExtra = vertcat(CMLExtra{:});
    end
%deals with ML
    function [MLSummary, MLTime, MLExtra] = MLAnalyze(ML, alignEvent)
        %Basic Analysis
        MLSummary = resultsIndTimeCourseMod(ML, false);
        %Time Analysis
        timeProcIndx = ML.getProcessIndex('TimePoints');
        alignIndx = ML.processes_{timeProcIndx}.getIndex(alignEvent);
        offSet = datenum(ML.processes_{timeProcIndx}.times_{alignIndx});
        MLTime = cellfun(@(x) (datenum(x.acquisitionDate_) - offSet) .* 1440, ML.movies_, 'UniformOutput', false);
        MLTime = vertcat(MLTime{:});
        %Conditional (Extra) Analysis
        MLExtra = cellfun(@MDAnalyze, ML.movies_, 'UniformOutput', false);
        MLExtra = vertcat(MLExtra{:});
    end
%deals with MD
    function [MDExtra] = MDAnalyze(MD)
        %Need blank if not used
        MDExtra.blank = [];
        %loads 'partitionResult'
        if analysisPara.doPartition
            load(MD.processes_{channel, MD.getProcessIndex('PartitionAnalysisProcess')}.outFilePaths_{1});
        end
        %loads 'tracks' and 'diffAnalysisRes'
        if analysisPara.doPartition
            load(MD.processes_{channel, MD.getProcessIndex('MotionAnalysisProcess')}.outFilePaths_{1});
        end
        %--------------------TrackPartitioning analysis---------------
        if analysisPara.doPartition
            %Column 1 : immobile
            %Column 2 : confined
            %Column 3 : free
            %Column 4 : directional
            %Column 5 : undetermined
            %Column 6 : determined
            %Column 7 : total
            %Initialize
            partitionSum = [0 0 0 0 0];
            trackTotal = [0 0 0 0 0];
            arrayfun(@(x, y) partitionTrack(x.classification(:,2), y.partitionFrac), tracks, partitionResult);
            MDExtra.partitionFrac = partitionSum ./ trackTotal;
            MDExtra.partitionFrac(6) = sum(partitionSum(1:4)) ./ sum(trackTotal(1:4));
            MDExtra.partitionFrac(7) = sum(partitionSum) ./ sum(trackTotal);
        end
        %Nested function that looks at each track
        function partitionTrack(diffClass, partition)
            nSubTrack = numel(partition);
            for iSubTrack = 1:nSubTrack
                if diffClass(iSubTrack) == 0
                    trackTotal(1) = trackTotal(1) + 1;
                    partitionSum(1) = partitionSum(1) + partition(iSubTrack);
                elseif diffClass(iSubTrack) == 1
                    trackTotal(2) = trackTotal(2) + 1;
                    partitionSum(2) = partitionSum(2) + partition(iSubTrack);
                elseif diffClass(iSubTrack) == 2
                    trackTotal(3) = trackTotal(3) + 1;
                    partitionSum(3) = partitionSum(3) + partition(iSubTrack);
                elseif diffClass(iSubTrack) == 3
                    trackTotal(4) = trackTotal(4) + 1;
                    partitionSum(4) = partitionSum(4) + partition(iSubTrack);
                elseif isnan(diffClass(iSubTrack))
                    trackTotal(5) = trackTotal(5) + 1;
                    partitionSum(5) = partitionSum(5) + partition(iSubTrack);
                end
            end
        end
        %---------------------
    end
%% Format Change
%Using resultsCombTimeCourseMod.m to store data in correct format
summary = cellfun(@resultsCombTimeCourseMod, summary, 'UniformOutput', false);
%adds time and name (because that's not in resultsCombTimeCourseMod.m)
for iCML = 1:nCML
    summary{iCML}.time = time{iCML};
    summary{iCML}.name = CMLs(iCML).name_;
end
%adds extra analysis if applicable
if analysisPara.doPartition
    for iCML = 1:nCML
        summary{iCML}.partitionFrac = extra{iCML}.partitionFrac;
    end
end
%% Sort
%order data from earliest time point to latest
for iCML = 1:nCML
    [summary{iCML}.time, sortIndx] = sort(summary{iCML}.time);
    summary{iCML}.numAbsClass = summary{iCML}.numAbsClass(sortIndx, :);
    summary{iCML}.numNorm0Class = summary{iCML}.numNorm0Class(sortIndx, :);
    summary{iCML}.probClass = summary{iCML}.probClass(sortIndx, :);
    summary{iCML}.diffCoefClass = summary{iCML}.diffCoefClass(sortIndx, :);
    summary{iCML}.confRadClass = summary{iCML}.confRadClass(sortIndx, :);
    summary{iCML}.ampClass = summary{iCML}.ampClass(sortIndx, :);
    summary{iCML}.ampNormClass = summary{iCML}.ampNormClass(sortIndx, :);
    summary{iCML}.ampStatsF20 = summary{iCML}.ampStatsF20(sortIndx, :);
    summary{iCML}.ampStatsL20 = summary{iCML}.ampStatsL20(sortIndx, :);
    summary{iCML}.rateMS = summary{iCML}.rateMS(sortIndx, :);
    if analysisPara.doPartition
        summary{iCML}.partitionFrac = summary{iCML}.partitionFrac(sortIndx, :);
    end
end
%% Plot and Save
timeCourseAnalysis_StandAlone(summary, outputDir, 'smoothingPara', smoothingPara, 'showPartitionAnalysis', analysisPara.doPartition);
save([outputDir filesep 'analysisData.mat'], 'directory_CML', 'analysisPara', 'summary');
end

