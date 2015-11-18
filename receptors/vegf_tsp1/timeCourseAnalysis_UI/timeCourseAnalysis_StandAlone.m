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
%       .getTimes           : function handle that returns commonInfo.times
%                             Don't rely on this!!!
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
ip.addParameter('showPartitionAnalysis', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('smoothingPara', .1, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParameter('nBootstrp', 100, @isnumeric);
ip.addParameter('timeResolution', 1, @isnumeric);
ip.addParameter('curveCompareAlpha', 0.05, @(x) isnumeric(x) && x>0 && x<1);
ip.addParameter('compareCurves', true, @(x) islogical(x)||isnumeric(x));
ip.addParameter('shiftPlotPositive', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('shiftTime', [], @(x) isnumeric(x));
ip.parse(varargin{:});
params = struct('showPartition', ip.Results.showPartitionAnalysis, 'smoothingPara', ip.Results.smoothingPara, 'nBootstrp', round(ip.Results.nBootstrp),...
    'timeResolution', ip.Results.timeResolution, 'curveCompareAlpha', ip.Results.curveCompareAlpha, 'compareCurves', ip.Results.compareCurves,...
    'shiftPlotPositive', ip.Results.shiftPlotPositive, 'shiftTime', ip.Results.shiftTime);
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
nConditions = numel(data);
timeShift = zeros(1, nConditions);
times = cellfun(@(x) x.time, data, 'UniformOutput', false);
names = cellfun(@(x) x.name, data, 'UniformOutput', false);
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
condColorAll = {...
    [0 0 0],... %black 1
    [0 0 1],... %blue 2
    [0 0.5 0],... %dark green 3
    [1 0 0],... %red 4
    [1 0.6 0.3],... %orange 5
    [0 1 0],... %green 6
    [1 0 1],... %magenta 7
    [0 1 1],... %cyan 8
    [1 1 0],... %yellow 9
    [0.7 0.5 0],... %brown 10
    [0.7 0.7 0.7]... %gray 11
    [0.5 0.5 1],... %purple 12
    [0.3 0.8 1],... %light blue 13
    };
colors = condColorAll(mod(1:nConditions, 13) + 1);
colors = reshape(colors,size(data));
%For saving plot data
figureData(37) = struct('titleBase', [], 'titleVariable', [], 'figureDir', [], 'fitData', [], 'data', [], 'getTimes', [], 'fitError', [], 'fitCompare', []);
commonInfo = struct('times', [], 'compareTimes', [], 'conditions', [], 'parameters', [], 'fullPath', [outputDir filesep 'figureData.mat'], 'timeShift', timeShift);
commonInfo.times = times;
commonInfo.conditions = names;
commonInfo.parameters = params;
iFD = 0;

%% Plot and Fit Smoothing Spline
%Initialize
defCond = {'immobile', 'confined', 'free', 'directed', 'undetermined', 'determined', 'total'};
%progressText
fprintf('Plotting figures: scatter plots\n');
%Each line calls the nested function plotData
%the first input subData must be cellarray of arrays
%So using cell fun convert data which is a cellarray of structure of arrays
%into a cellarray of arrays
%Other inputs are explained in plotData
plotFigure(cellfun(@(x) x.numAbsClass, data, 'UniformOutput', false), 'Absolute Number of Class Types', defCond, 'Number of tracks (molecules)', true);
plotFigure(cellfun(@(x) x.numNorm0Class, data, 'UniformOutput', false), 'Normalized Number of Class Types', defCond, 'Relative number of tracks', true);
plotFigure(cellfun(@(x) x.probClass, data, 'UniformOutput', false), 'Probability of Class Types', {'immobile', 'confined', 'free', 'directed', 'determined'}, 'Probability', true);
plotFigure(cellfun(@(x) x.diffCoefClass, data, 'UniformOutput', false), 'Diffusion Coefficient', defCond(1:4), 'Diffusion coefficient(pixels^2/frame)', true); %no 5th
plotFigure(cellfun(@(x) x.confRadClass, data, 'UniformOutput', false), 'Confinement Radius', defCond(1:2), 'Radius (pixels)', true);%no 3 4 5th column
plotFigure(cellfun(@(x) x.ampClass, data, 'UniformOutput', false), 'Fluorescence Amplitude', defCond(1:5), 'Intensity (arbitrary units)', true);
plotFigure(cellfun(@(x) x.ampNormClass, data, 'UniformOutput', false), 'Normalized Fluorescence Amplitude', defCond(1:5), 'Normalized intensity (monomer units)', true);
%plotData(cellfun(@(x) x.ampStatsF20, data, 'UniformOutput', false), '', {'mean', 'first mode mean', 'first mode std', 'first mode fraction', 'number of modes', 'normalized mean'}, '', true);
%plotData(cellfun(@(x) x.ampStatsL20, data, 'UniformOutput', false), '', {'mean', 'first mode mean', 'first mode std', 'first mode fraction', 'number of modes', 'normalized mean'}, '', true);
plotFigure(cellfun(@(x) x.rateMS, data, 'UniformOutput', false), 'Merging and Spliting', {'merging', 'spliting'}, '(per frame)', true);
%Do only if input specify that this plot be shown. Will cause error if
%data.partitionFrac is not present
if params.showPartition
    plotFigure(cellfun(@(x) x.chemEnergy, data, 'UniformOutput', false), 'Chemical Energy of Localization', defCond(1:5), 'Chemical Energy (arbitrary energy units)', false);
    plotFigure(cellfun(@(x) x.locFreq, data, 'UniformOutput', false), 'Localization Frequency', defCond(1:5), 'k on (arbitrary units)', true);
    plotFigure(cellfun(@(x) x.delocFreq, data, 'UniformOutput', false), 'Delocalization Frequency', defCond(1:5), 'k off (arbitrary units)', true);
    plotFigure(cellfun(@(x) x.eqCond, data, 'UniformOutput', false), 'Equilibrium Condition', defCond(1:5), 'Proximity to equilibrium condition (arbitrary units)', false);
end
%get rid of figure data that was not plotted
mask = arrayfun(@(x) ~isempty(x.fitData), figureData);
figureData = figureData(mask);
pause(1);
%progressText
fprintf('\b Complete\n');

%% Nested function for plotting
% Splits data structure elements by columns and calls the plotting function
%In other words, converts subData which is cell array of arrays into cell
%array of columns.
%subData            : contains all data to be plotted
%title_Base         : part of plot title that does not change
%title_variable     : part of plot title that does change depending on the
%                     column number
%yLabelName         : used for ylabel of the plot
    function plotFigure(subData, title_Base, title_Variable, yLabelName, isYMin0)
        %initialization
        nColumns = numel(title_Variable);
        plotData(nColumns) = struct('fitData', [], 'condition', []);
        %determine maximum y value to determine y axis limit
        maxValue = max(cellfun(@(x) max(max(x(~(isnan(x)|isinf(x))))), subData));
        axisTick = 10^(round(log10(maxValue) + 0.2)-1);
        if isYMin0
            yMin = 0;
        else
            minValue = min(cellfun(@(x) min(min(x(~(isnan(x)|isinf(x))))), subData));
            axisTick = max(10^(round(log10(abs(minValue)) + 0.2)-1), axisTick);
            yMin = (ceil(minValue / axisTick) - 1) * axisTick;
        end
        yMax = (floor(maxValue / axisTick) + 1) * axisTick;
        %plot by column
        for iColumns = 1:nColumns
            subSubData = cellfun(@(x) x(:,iColumns), subData, 'UniformOutput', false);
            plotTitle = [title_Base ' ' title_Variable{iColumns}];
            [fitData] = plotMultipleSmoothingSpline(outputDirFig, subSubData, times, names, colors, plotTitle, yLabelName, params.smoothingPara, yMax, yMin, timeShift);
            plotData(iColumns).fitData = fitData;
            plotData(iColumns).condition = title_Variable{iColumns};
            %Save data in figureData
            iFD = iFD + 1;
            figureData(iFD).titleBase = title_Base;
            figureData(iFD).titleVariable = title_Variable{iColumns};
            figureData(iFD).fitData = fitData;
            figureData(iFD).data = subSubData;
            figureData(iFD).figureDir = [outputDirFig filesep plotTitle '.fig'];
            figureData(iFD).getTimes = @() commonInfo.times;
            figureData(iFD).yMax = yMax;
            figureData(iFD).yMin = yMin;
            figureData(iFD).yLabel = yLabelName;
        end
    end

%% BootStrap Analysis of Fits
%This is used to determine the standard error of the fitted curve
%determination of time limit for the analysis--------------------------
%determine the limit of each conditions

if(~params.nBootstrp)
    % If nBootstrp is zero, then assign empty values and quit
    commonInfo.analysisTimes = [];
    if(~isempty(figureData))
        [figureData.fitError] = deal([]);
    end
    return;
end

timeMax = cellfun(@(x) max(x), commonInfo.times);
timeMin = cellfun(@(x) min(x), commonInfo.times);
%determine the overall range of all conditions (This is useful later when comparing two curves)
timeMaxMax = max(timeMax);
timeMinMin = min(timeMin);
%convert to number divisible by timeResolution
timeMinMin = timeMinMin - mod(timeMinMin, params.timeResolution);
timeMaxMax = timeMaxMax - mod(timeMaxMax, -params.timeResolution);
analysisTime_Union = timeMinMin : params.timeResolution : timeMaxMax;
%Determine the indx and value of timemin and max for each conditions
timeLimit = cell(1, nConditions);
timeLimitIndx = cell(1, nConditions);
for iCond = 1:nConditions
    timeLimitIndx{iCond} = [find(analysisTime_Union <= timeMin(iCond), 1, 'last'), find(analysisTime_Union >= timeMax(iCond), 1, 'first')];
    timeLimit{iCond} = [analysisTime_Union(timeLimitIndx{iCond}(1)), analysisTime_Union(timeLimitIndx{iCond}(2))];
end
timeLimit = reshape(timeLimit,size(data));
%store this in commonInfo
commonInfo.analysisTimes = cellfun(@(x) x(1) : params.timeResolution : x(2), timeLimit, 'UniformOutput', false);
%for progress display
nFig = numel(figureData);
progressTextMultiple('Determining confidence interval', nFig);
%call determineSE_Bootstrp.m
fitError = arrayfun(@(x) determineSE(x.data, commonInfo.times, params.nBootstrp, params.timeResolution, timeLimit, params.smoothingPara), figureData, 'Uniformoutput', false, 'ErrorHandler', @determineSEEH);
[figureData.fitError] = fitError{:};

%% Add Standard Error to Figures
timeCourseAnalysis.plot.standardErrorFigure(commonInfo,figureData,true,outputDirFig2);
% lineObj = cell(1,nConditions);
% for iFig = 1:nFig
%     %open up the figure and make it current figure
%     %figObj = openfig(figureData(iFig).figureDir);
%     figObj = figure();
%     hold on;
%     %plot standard errors
%     for iCond = 1:nConditions
%         fitValues = figureData(iFig).fitData{iCond}(commonInfo.analysisTimes{iCond});
%         fitValues = fitValues';
%         plot(commonInfo.analysisTimes{iCond}, fitValues + figureData(iFig).fitError{iCond}, ':', 'color', colors{iCond});
%         plot(commonInfo.analysisTimes{iCond}, fitValues - figureData(iFig).fitError{iCond}, ':', 'color', colors{iCond});
%         lineObj{iCond} = plot(commonInfo.analysisTimes{iCond}, fitValues, 'color', colors{iCond});
%         plot([timeShift(iCond), timeShift(iCond)], [figureData(iFig).yMax, figureData(iFig).yMin], 'Color', colors{iCond});
%     end
%     lineObj2 = [lineObj{:}];
%     legend(lineObj2, commonInfo.conditions);
%     title([figureData(iFig).titleBase, ' ', figureData(iFig).titleVariable])
%     xlabel('Time (min)');
%     ylabel(figureData(iFig).yLabel);
%     ylim([figureData(iFig).yMin, figureData(iFig).yMax]);
%     %save and close
%     savefig(figObj, [outputDirFig2, filesep, figureData(iFig).titleBase, ' ', figureData(iFig).titleVariable, '.fig']);
%     close(figObj);
% end
pause(1);

%% Compare Fitted Curves
%progress display
progressTextMultiple('Comparing fittted curves', nFig);
%pairs up curve indices for comparison----------------
%And
%determines the commonInfo.compareTimes
%which are simply intersection of two commonInfo.analysisTimes sets
%And
%determines the indx of these overlap------------
nPair = nConditions * (nConditions - 1) / 2;
curveIndxPair = cell(1, nPair);
commonIndx1 = cell(1, nPair);
commonIndx2 = cell(1, nPair);
iPair_ = 0;
for indx1 = 1:nConditions-1
    for indx2 = indx1+1:nConditions
        iPair_ = iPair_ + 1;
        curveIndxPair{iPair_} = [indx1, indx2];
        commonMinIndx = max(timeLimitIndx{indx1}(1), timeLimitIndx{indx2}(1));
        commonMaxIndx = min(timeLimitIndx{indx1}(2), timeLimitIndx{indx2}(2));
        %for indx1
        commonIndx1{iPair_} = [commonMinIndx, commonMaxIndx] - timeLimitIndx{indx1}(1) + 1;
        %for indx2
        commonIndx2{iPair_} = [commonMinIndx, commonMaxIndx] - timeLimitIndx{indx2}(1) + 1;
        %indx to value conversion
        commonInfo.compareTimes{iPair_} = commonInfo.analysisTimes{indx1}(commonIndx1{iPair_}(1):commonIndx1{iPair_}(2));
    end
end
%determine p-values at each timepoints = commonInfo.compareTimes
%deals out figureData to callGetP
fitCompare = arrayfun(@callGetP, figureData, 'UniformOutput', false);
[figureData.fitCompare] = fitCompare{:};

%% Nested Function: deals out pair of curves to getP
    function [fitCompare]  = callGetP(FD)
        if isempty(curveIndxPair)
            pValue = [];
        else
            pValue = cellfun(@(x, y, index1, index2) getPValues(FD.fitData{x(1)}, FD.fitData{x(2)}, FD.fitError{x(1)}(index1(1):index1(2)), FD.fitError{x(2)}(index2(1):index2(2)), y), curveIndxPair, commonInfo.compareTimes, commonIndx1, commonIndx2, 'UniformOutput', false, 'ErrorHandler', @getPEH);
        end
        fitCompare(nConditions, nConditions) = struct('geoMeanP', [], 'p', [], 'timeIndx', []);
        for iPair = 1:nPair
            fitCompare(curveIndxPair{iPair}(1), curveIndxPair{iPair}(2)).p = pValue{iPair};
            fitCompare(curveIndxPair{iPair}(1), curveIndxPair{iPair}(2)).geoMeanP = geomean(pValue{iPair});
            fitCompare(curveIndxPair{iPair}(1), curveIndxPair{iPair}(2)).timeIndx = iPair;
        end
        progressTextMultiple();
    end

%% Save
save(commonInfo.fullPath, 'commonInfo', 'figureData');

end
%% Local Functions
%function for progressDisplay and calls determineSE_Bootstrp.m
function [fitError] = determineSE(data, time, nBoot, timeResolution, timeLimit, smoothingPara)
fitError = determineSmoothSplineSE(data, time, nBoot, timeResolution, timeLimit, smoothingPara);
progressTextMultiple();
end
%error handle for determineSE_Bootstrp.m
function [fitError] = determineSEEH(varargin)
fitError = [];
warning(['Standard error determination for figure ' num2str(varargin{1}.index) ' has failed']);
error(varargin{1});
end
%determines p-values
function [pValue] = getPValues(fitData1, fitData2, fitError1, fitError2, times)
pValue = arrayfun(@(x, y, t) getP(fitData1(t), fitData2(t), x, y), fitError1, fitError2, times);
end
%determines p-value using z-test
function [pValue] = getP(mean1, mean2, SE1, SE2)
z = abs(mean1 - mean2) ./ sqrt(SE1.^2 + SE2.^2);
pValue = 1 - diff(normcdf([-z,z]));
end
%error handle for p value determination
function [] = getPEH(varargin)
warning(['Curve Comparison for figure ' num2str(varargin{1}.index) ' has failed']);
error(varargin{1});
end
