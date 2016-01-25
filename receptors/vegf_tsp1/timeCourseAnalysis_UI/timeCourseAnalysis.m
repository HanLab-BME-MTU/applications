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
%       'channels'               : Channel of MD to be analyzed. Default = 1
%       'doNewAnalysis'           : (logical) true: always do new analysis even if
%                                 the analysis has already been done.
%                                 false: avoid doing the analysis again if
%                                 analysis has already been done and the
%                                 analysis parameters are identical.
%                                 true by default
%       'doPartitionAnalysis'   : (logical) to analyzes trackPartitioning
%                                 process or not. false by default
%       'shiftPlotPositive'     : (logical) determines whether or not whole
%                                 plots are shifted so that no negative
%                                 time values are present. In other words,
%                                 minimum time point is taken to be the
%                                 zero. Default is false. <don't use>
%       'start2zero'            : (logical) sets the average start time to
%                                 be zero, even if the align event is not
%                                 'start'
%       'channelNames'          : (cellstr) cell array of channel names
%        detectOutliers_k_sigma : (numeric, scalar) see detectOutliers
%
%
%Tae H Kim, July 2015

%% Input
%convert CMLs if it's cell of strings
nCML = numel(CMLs);
if iscellstr(CMLs)
    directory_CML = CMLs;
    clear CMLs
    for iCML = 1:nCML
        fprintf('Loading CombinedMovieList %g/%g\n', iCML, nCML);
        CMLs(iCML) = CombinedMovieList.load(directory_CML{iCML});
    end
elseif isa(CMLs, 'CombinedMovieList')
    directory_CML = arrayfun(@(x) x.getFullPath(), CMLs, 'UniformOutput', false); %#ok<NASGU>
else
    error('CMLs must be class object CombinedMovieList');
end
%input parser
ip = inputParser;
ip.CaseSensitive = false;
% Extra (Unmatched) parameters forwarded to timeCourseAnalysis_StandAlone
ip.KeepUnmatched = true;
ip.StructExpand = true;
ip.addRequired('outputDir', @ischar);
ip.addParameter('channels', [], @isnumeric);
ip.addParameter('doPartition', false, @(x) isnumeric(x) || islogical(x));
ip.addParameter('doNewAnalysis', true, @(x) isnumeric(x) || islogical(x));
ip.addParameter('start2zero', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('channelNames', false, @(x) iscellstr(x));
ip.parse(outputDir, varargin{:});
%for easier access to ip.Result variables
analysisPara = ip.Results;
%checks if CMLs are loaded or not
%AND determines how many ML there are total
nMLTot = 0;
for iCML = 1:nCML
    progressDisp = fprintf('Loading Combined Movie List %g/%g\n', iCML, nCML); %loading progress display
    % if not loads it
    if numel(CMLs(iCML).movieLists_) ~= numel(CMLs(iCML).movieListDirectory_)
        CMLs(iCML).sanityCheck();
    end
    %checks if all MD has necessary processes
    arrayfun(@(x) timeCourseAnalysis.MLCheck(x, analysisPara), CMLs(iCML).movieLists_);
    fprintf(repmat('\b', 1, progressDisp)); %loading progress display
    nMLTot = nMLTot + numel(CMLs(iCML).movieLists_);
end

%% Main Time Course Analysis (CML-level)
% Analyze all MovieData in parallel first
timeCourseAnalysis.analyzeMDsInParallel(CMLs,analysisPara.doNewAnalysis);
%Progress Display
progressTextMultiple('Time Course Analysis', nMLTot);
%Using resultsIndTimeCourseMod.m to do basic analysis
%and extract time data and align
[summary, time, extra, startTime] = arrayfun(@(CML) timeCourseAnalysis.CMLAnalyze(CML,analysisPara), CMLs, 'UniformOutput', false);
startTime = [startTime{:}];

%% Format Change
% Break summary apart by channel
summary = cellfun(@(s) num2cell(s,1),summary,'UniformOutput',false);
% Should a nCML x nChannel cell array
% Made an assumption that all that CMLs have the same channel length here
summary = vertcat(summary{:});
%Using resultsCombTimeCourseMod.m to store data in correct format
summary = cellfun(@resultsCombTimeCourseMod, summary, 'UniformOutput', false);
%adds time and name (because that's not in resultsCombTimeCourseMod.m)
if(isempty(analysisPara.channels))
    analysisPara.channels = 1:size(summary,2);
end
for iCML = 1:nCML
    for iChannel = analysisPara.channels;
        summary{iCML,iChannel}.time = time{iCML};
        summary{iCML,iChannel}.name = analysisPara.channelNames{iChannel};
        if(~isempty(CMLs(iCML).name_))
            if(~isempty(summary{iCML,iChannel}.name))
                summary{iCML,iChannel}.name = [ '/' summary{iCML,iChannel}.name];
            end
            summary{iCML,iChannel}.name = [ CMLs(iCML).name_ summary{iCML,iChannel}.name];
        end
    end
end
% HACK, FIX
% essentially let the rest of the analysis use linear indexing to access
% summary
nCML = numel(summary);
%adds extra analysis if applicable
if analysisPara.doPartition
    for iCML = 1:nCML
        summary{iCML}.chemEnergy = vertcat(extra{iCML}.chemEnergy);
        summary{iCML}.locFreq = vertcat(extra{iCML}.locFreq);
        summary{iCML}.delocFreq = vertcat(extra{iCML}.delocFreq);
        summary{iCML}.eqCond = vertcat(extra{iCML}.eqCond);
    end
end
%% Sort
%order data from earliest time point to latest
for iCML = 1:nCML
    [~, sortIndx] = sort(summary{iCML}.time);
    % Name should not be sorted
    name = summary{iCML}.name;
    summary{iCML} = structfun(@(x) x(sortIndx,:),summary{iCML},'UniformOutput',false,'ErrorHandler',@(~,x) x);
    summary{iCML}.name = name;
end

%% Shift time
shiftTime = [];
if analysisPara.start2zero
    shiftTimeIndx = startTime~=0;
    offset = - mean(startTime(shiftTimeIndx));
    shiftTime = zeros(1, nCML);
    for iCML = find(shiftTimeIndx)
        shiftTime(iCML) = offset;
    end
end

%% Plot by Calling StandAlone Function
timeCourseAnalysis_StandAlone(summary, outputDir, ...
      ip.Unmatched ...
    , 'shiftTime', shiftTime ...
    );
%% Save
save([outputDir filesep 'analysisData.mat'], 'directory_CML', 'analysisPara', 'summary');
end

