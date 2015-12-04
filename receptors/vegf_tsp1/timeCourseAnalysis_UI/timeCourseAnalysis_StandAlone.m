function [commonInfo, figureData] = timeCourseAnalysis_StandAlone(data, outputDir, varargin)
%Time course analysis of Movie Data
%
%SYNOPSIS function [] = TimeCourseAnalysis_StandAlone(data, outputDir, varargin)
%
%INPUT
%    data           : Cell of array of structure with timecourse analyzed MD
%       .name
%       .numAbsClass    : Absolute number of particles in the various
%                         motion classes.
%                         Row = time points in time course (see time).
%                         Columns = immobile, confined, free, directed,
%                         undetermined, determined, total.
%       .numNorm0Class  : Normalized number of particles in the various
%                         motion classes, such that the first time = 1.
%                         Rows and columns as numAbsClass.
%       .probClass      : Probability of the various motion classes.
%                         Rows as above.
%                         Columns = immobile, confined, free, directed (all
%                         relative to determined); and determined relative
%                         to total.
%       .diffCoefClass  : Mean diffusion coefficient in the various
%                         motion classes.
%                         Rows as above.
%                         Columns = immobile, confined, free, directed,
%                         undetermined.
%       .confRadClass   : Mean confinement radius in the various
%                         motion classes.
%                         Rows and columns as above.
%       .ampClass       : Mean amplitude of particles in the various motion
%                         classes.
%                         Rows and columns as above.
%       .ampNormClass   : Mean normalized amplitude of particles in the
%                         various motion classes.
%                         Rows and columns as above.
%       .ampStatsF20    : Amplitude statistics for particles in first 20
%                         frames.
%                         Rows as above.
%                         Columns = mean, first mode mean, first mode std,
%                         first mode fraction, number of modes, normalized
%                         mean.
%       .ampStatsL20    : Amplitude statistics for particles in last 20
%                         frames.
%                         Rows and columns as above.
%       .rateMS         : Rate of merging and rate of splitting per
%                         feature. Columns: merging, splitting.
%       .time           : List of relative time points (single column)
% 
%       *following elements may or may not be present
% 
%       .partitionFrac  : Probability of finding the various motion classes
%                         within a mask
%                         Row = time points in time course (see time).
%                         Columns = immobile, confined, free, directed,
%                         undetermined, determined, total.
% 
%   outputDir       : file location where all figures and data will be saved
% 
%   varargin        : name_value pair
%
%       showPartitionAnalysis   : logical determining if .partitionFrac will
%                                be shown or not. Can be 0 or 1. If true,
%                                must have data.partitionFrac
%       smoothingPara           : parameter used for smoothing spline fit
%       nBootstrp               : number of bootstrap data sets to use for
%                                 bootstrap analysis for standard error
%                                 determination. default = 100
%       timeResolution          : time resolution of bootstrap analysis
%                                 (in minutes). Default is 1;
%       curveCompareAlpha       : cut off p value needed to reject the null
%                                 hypothesis that two curves are the same.
%                                 default = 0.05
%       compareCurves           : (logical) determines if the program
%                                 should compare curves within a figure to
%                                 see if they are significantly different
%                                 or not.
%       detectOutliers_k_sigma  : (numeric, scalar) see detectOutliers
% 
%OUTPUT
%   commonInfo      : contains information and data common to all figures
%                     or plots
%       .times          : Cellarray of relative time points in a single
%                         column that correspond to data points.
%                         The elements correspond to 'commonInfo.conditions'.
%       .analysisTimes  : Cellarray of time columns where fit was analyzed
%                         for SE (or 1 sigma confidence interval)
%       .compareTimes   : Cellarray of time columns where two fit were
%                         compared for similarity
%       .conditions     : name or the conditions plotted (ie. VEGF- AAL-)
%       .params         : parameters used for analysis
%           .showPartition      : logical determining if .partitionFrac will
%                                be shown or not. Can be 0 or 1. If true,
%                                must have data.partitionFrac
%           .smoothingPara      : parameter used for smoothing spline fit
%           .nBootstrp          : number of bootstrap data sets to use for
%                                 bootstrap analysis for standard error
%                                 determination. default = 100
%           .timeResolution     : time resolution of bootstrap analysis
%                                 (in minutes). Default is 1;
%           .curveCompareAlpha  : cut off p value needed to reject the null
%                                 hypothesis that two curves are the same.
%                                 default = 0.05
%           .compareCurves      : (logical) determines if the program
%                                 should compare curves within a figure to
%                                 see if they are significantly different
%                                 or not. Default is true. -Unused-
%           .shiftPlotPositive  : (logical) determines whether or not whole
%                                 plots are shifted so that no negative
%                                 time values are present. In other words,
%                                 minimum time point is taken to be the
%                                 zero. Default is false.
%           .shiftTime          : (array, numeric) how much the aligned time are
%                                 shifted. if ==1, align Event will be at 1
%       .fullPath       : fullpath of where commonInfo and figureData are
%                         saved
%       .timeShift      : how much time values have been changed from the
%                         original values due to time shifting mechanisms
%                         like .start2zero
%                         To get to the original value, simply subtract
%                         .timeShift from the recorded time values
%                         Here, original value means that the align event is
%                         at 0.
%   figureData      : array of information and data specific to a figure
%       .titleBase          : constant part of plot title
%       .titleVariable      : variable part of plot title
%       .figureDir          : full path of related figure
%
%       *Following structure elements contain a cellarray or an array.
%       *The elements of these cellarrays or arrays correspond to 'commonInfo.conditions'.
%
%       .fitData            : cellarray of smoothingSpline fit
%       .data               : cellarray of plotted data sets
%                             (each element in a column)
%       .yMax               : y-axis maximum on plot
%       .yMin               : y-axis minimum on plot
%       .yLabel             : y axis label
%       .fitError           : cellarray of the standard error calculated
%                             from bootstrap analysis
%       .fitCompare(i,j)    : n by n array that compares the curve
%                             fits to each other. Not empty only if i<j.
%                             If i>=j, the strucutre elements will be empty.
%           .geoMeanP       : geometric mean of p-values
%           .p              : list of all p-values. Higher p-value means
%                             null hypothesis is more likely.
%           .timeIndx       : the index of commonInfo.compareTimes that
%                             correspond to this fit compare
%                             
%
%Tae Kim, July 2015

%% Input
%input parse
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.StructExpand = true;
ip.addParameter('showPartitionAnalysis', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('smoothingPara', .01, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParameter('nBootstrp', 100, @isnumeric);
ip.addParameter('timeResolution', 1, @isnumeric);
ip.addParameter('curveCompareAlpha', 0.05, @(x) isnumeric(x) && x>0 && x<1);
ip.addParameter('compareCurves', true, @(x) islogical(x)||isnumeric(x));
ip.addParameter('shiftPlotPositive', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('shiftTime', [], @(x) isnumeric(x));
ip.addParameter('detectOutliers_k_sigma', [], @(x) isnumeric(x));
ip.parse(varargin{:});
params = ip.Results;
params.showPartition = params.showPartitionAnalysis;

outputDirFig = [outputDir filesep 'figures'];
outputDirFig2 = [outputDir filesep 'figures_SE'];
%makes sure outputDir folder exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(outputDirFig, 'dir')
    mkdir(outputDirFig);
end
if ~exist(outputDirFig2, 'dir')
    mkdir(outputDirFig2);
end

%% Inititalization
data = cell2mat(data);
nConditions = numel(data);
timeShift = zeros(1, nConditions);
% extractField = @(field) cellfun(data.(field),'UniformOutput',false);
times ={data.time}';
names = {data.name}';
%Shift values
if params.shiftPlotPositive
    minTime = min(cellfun(@min, times));
    if minTime < 0
        times = cellfun(@(x) x - minTime, times, 'UniformOutput', false);
        timeShift = timeShift - minTime;
    end
end
if numel(params.shiftTime) == nConditions
    for iCond = 1:nConditions
        times{iCond} = times{iCond} + params.shiftTime(iCond);
        timeShift(iCond) = timeShift(iCond) + params.shiftTime(iCond);
    end
end
%Color initialization
% colors = timeCourseAnalysis.plot.getColors(data);
%For saving plot data
% figureData(37) = struct('titleBase', [], 'titleVariable', [], 'figureDir', [], 'fitData', [], 'data', [], 'getTimes', [], 'fitError', [], 'fitCompare', []);
figureData = {};
commonInfo = struct('times', [], 'compareTimes', [], 'conditions', [], 'parameters', [], 'fullPath', [outputDir filesep 'figureData.mat'], 'timeShift', timeShift);
commonInfo.times = times;
commonInfo.conditions = names;
commonInfo.parameters = params;
commonInfo.outputDirFig = outputDirFig;
commonInfo.outputDirFig2 = outputDirFig2;
% iFD = 0;

%% Plot and Fit Smoothing Spline
%Initialize
defCond = { ...
      'immobile' ...      % 1
    , 'confined' ...      % 2 
    , 'free' ...          % 3
    , 'directed' ...      % 4
    , 'undetermined' ...  % 5
    , 'determined' ...    % 6
    , 'total' ...         % 7
    , 'sub-diffusive' ... % 8
    };
ampLabels = { ...
    'Fluorescence Amplitude Overall (a.u.)' ... % 1
    ,'First Mode Mean (a.u.)', 'First Mode Std (a.u.)' ... % 2
    , 'First Mode Fraction' ... % 3
    , 'Number of Modes' ... % 4
    ,'Normalized Fluorescence Amplitude Overall (monomer units)' ... % 5
    };

commonInfo.defCond = defCond;
commonInfo.ampLabels = ampLabels;

%progressText
fprintf('Plotting figures: scatter plots\n');

%Each line calls the nested function plotData
%the first input subData must be cellarray of arrays
%So using cell fun convert data which is a cellarray of structure of arrays
%into a cellarray of arrays
%Other inputs are explained in plotData
calcFigure({data.numAbsClass}', 'Absolute Number of Class Types', ...
    defCond, 'Number of tracks (molecules)');
calcFigure({data.numNorm0Class}', 'Normalized Number of Class Types', ...
    defCond, 'Number of tracks relative to time 0');
calcFigure({data.densityAbsClass}', 'Absolute Density of Class Types', ...
    defCond, 'Density of tracks (molecules/pixel^2)');
calcFigure({data.densityNorm0Class}', 'Normalized Density of Class Types', ...
    defCond, 'Density of tracks relative to time 0');
calcFigure({data.probClass}', 'Probability of Class Types', ...
    defCond([1:4 6 8]), ...
    'Probability');

calcFigure({data.diffCoefClass}', 'Diffusion Coefficient', ...
    defCond(1:4), 'Diffusion coefficient (pixels^2/frame)'); %no 5th
calcFigure({data.confRadClass}', 'Confinement Radius', ... 
    defCond(1:2), 'Confinement radius (pixels)');%no 3 4 5th column

calcFigure({data.ampClass}', 'Fluorescence Amplitude', ...
    defCond(1:5), 'Intensity (arbitrary units)');
calcFigure({data.ampNormClass}', 'Normalized Fluorescence Amplitude', ... 
    defCond(1:5), 'Normalized intensity (monomer units)');
calcFigure({data.ampStatsF20}', 'First 20 Frames - ', ...
    ampLabels, '');
calcFigure({data.ampStatsL20}', 'Last 20 Frames - ', ...
    ampLabels, '');

calcFigure({data.rateMS}', 'M & S Rate', {'merging', 'splitting'}, '(per frame per particle)');
calcFigure({data.msTimeInfo}', 'M & S Time Information', ...
    {'merge-to-split time', 'split-to-merge (self) time','split-to-merge (other) time','merge-to-end time','start-to-split time'}, '(frames)');

%Do only if input specify that this plot be shown. Will cause error if
%data.partitionFrac is not present
if params.showPartition
    calcFigure({data.chemEnergy}', 'Chemical Energy of Localization', ...
        defCond(1:5), 'Chemical Energy (arbitrary energy units)', false);
    calcFigure({data.locFreq}', ...
        'Localization Frequency', defCond(1:5), 'k on (arbitrary units)');
    calcFigure({data.delocFreq}', 'Delocalization Frequency', ...
        defCond(1:5), 'k off (arbitrary units)');
    calcFigure({data.eqCond}', 'Equilibrium Condition', ...
        defCond(1:5), 'Proximity to equilibrium condition (arbitrary units)', false);
end

%get rid of figure data that was not plotted
figureData = [figureData{:}];
mask = arrayfun(@(x) ~any(cellfun('isempty',x.fitData)), figureData);
figureData = figureData(mask);

%% Save
save(commonInfo.fullPath, 'commonInfo', 'figureData');

%% Plot
timeCourseAnalysis.plot.scatterFigure(commonInfo,figureData,commonInfo.outputDirFig);
pause(1);
%progressText
fprintf('\b Complete\n');

%% Nested function for plotting
% Splits data structure elements by columns and sotre information for
% plotting
%In other words, converts subData which is cell array of arrays into cell
%array of columns.
%subData            : contains all data to be plotted
%title_Base         : part of plot title that does not change
%title_variable     : part of plot title that does change depending on the
%                     column number
%yLabelName         : used for ylabel of the plot
    function calcFigure(subData, title_Base, title_Variable, yLabelName, isYMin0)
        if(nargin < 4)
            yLabelName = '';
        end
        if(nargin < 5)
            isYMin0 = true;
        end
        figureData{end+1} = timeCourseAnalysis.calcFigureData(commonInfo, subData, title_Base, title_Variable, yLabelName, isYMin0);
    end


%% BootStrap Analysis of Fits
%This is used to determine the standard error of the fitted curve
%determination of time limit for the analysis--------------------------
%determine the limit of each conditions

if(~params.nBootstrp)
    % If nBootstrp is zero, then assign empty values and quit
    commonInfo.analysisTimes = [];
    commonInfo.timeLimitIndx = [];
    if(~isempty(figureData))
        [figureData.fitError] = deal([]);
    end
else

[commonInfo.analysisTimes, timeLimit, commonInfo.timeLimitIndx] = timeCourseAnalysis.getAnalysisTimes(commonInfo.times,params.timeResolution);
    
% Computer standard error
determineSEInParallel = true;
nFig = numel(figureData);
if(~determineSEInParallel)
    %for progress display
    
    progressTextMultiple('Determining confidence interval', nFig);
else
    warning('off','parallel:lang:spmd:RemoteTransfer');
    disp('Determining confidence interval');
    disp('Parallel progress not available');
    dFigureData = distributed(figureData);
end
%call determineSE_Bootstrp.m
fitError = arrayfun(@(x) determineSE(x.data, commonInfo.times, params.nBootstrp, params.timeResolution, timeLimit, params.smoothingPara, x.inOutFlag), dFigureData, 'Uniformoutput', false, 'ErrorHandler', @determineSEEH);
if(determineSEInParallel)
    fitError = gather(fitError);
end
[figureData.fitError] = fitError{:};

%% Add Standard Error to Figures
timeCourseAnalysis.plot.standardErrorFigure(commonInfo,figureData,true,outputDirFig2);

pause(1);

%% Compare Fitted Curves
[fitCompare, commonInfo.compareTime] = timeCourseAnalysis.compareFittedCurves(commonInfo, figureData);
[figureData.fitCompare] = fitCompare{:};


end % end of if(~params.nBootstrp)

%% Save
save(commonInfo.fullPath, 'commonInfo', 'figureData');

end
%% Local Functions
%function for progressDisplay and calls determineSE_Bootstrp.m
function [fitError] = determineSE(data, time, nBoot, timeResolution, timeLimit, smoothingPara, inOutFlag)
    fitError = determineSmoothSplineSE(data, time, nBoot, timeResolution, timeLimit, smoothingPara, inOutFlag);
    if(numlabs == 1)
        progressTextMultiple();
    end
end
%error handle for determineSE_Bootstrp.m
function [fitError] = determineSEEH(varargin)
fitError = [];
warning(['Standard error determination for figure ' num2str(varargin{1}.index) ' has failed']);
error(varargin{1});
end
