function [] = timeCourseAnalysis_StandAlone(data, outputDir, varargin)
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
% 
%Tae Kim, July 2015

%% Input
%input parse
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('showPartitionAnalysis', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('smoothingPara', .95, @(x) isnumeric(x) && x>=0 && x<=1);
ip.parse(varargin{:});
showPartition = ip.Results.showPartitionAnalysis;
smoothingPara = ip.Results.smoothingPara;
outputDirFig = [outputDir filesep 'figures'];
%makes sure outputDir folder exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(outputDirFig, 'dir')
    mkdir(outputDirFig);
end
%% Inititalization
nConditions = numel(data);
times = cellfun(@(x) x.time, data, 'UniformOutput', false);
names = cellfun(@(x) x.name, data, 'UniformOutput', false);
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
%% Plot
%Each line calls the nested function plotData
%the first input subData must be cellarray of arrays
%So using cell fun convert data which is a cellarray of structure of arrays
%into a cellarray of arrays
%Other inputs are explained in plotData
plotData(cellfun(@(x) x.numAbsClass, data, 'UniformOutput', false), 'Absolute Number of Class Types', {'immobile', 'confined', 'free', 'directed', 'undetermined', 'determined', 'total'}, 'Number of tracks (molecules)');
plotData(cellfun(@(x) x.numNorm0Class, data, 'UniformOutput', false), 'Normalized Number of Class Types', {'immobile', 'confined', 'free', 'directed', 'undetermined', 'determined', 'total'}, 'Relative number of tracks');
plotData(cellfun(@(x) x.probClass, data, 'UniformOutput', false), 'Probability of Class Types', {'immobile', 'confined', 'free', 'directed', 'determined'}, 'Probability');
plotData(cellfun(@(x) x.diffCoefClass, data, 'UniformOutput', false), 'Diffusion Coefficient', {'immobile', 'confined', 'free', 'directed'}, 'Diffusion coefficient(pixels^2/frame)'); %no 5th
plotData(cellfun(@(x) x.confRadClass, data, 'UniformOutput', false), 'Confinement Radius', {'immobile', 'confined'}, 'Radius (pixels)');%no 3 4 5th column
plotData(cellfun(@(x) x.ampClass, data, 'UniformOutput', false), 'Fluorescence Amplitude', {'immobile', 'confined', 'free', 'directed', 'undetermined'}, 'Intensity (arbitrary units)');
plotData(cellfun(@(x) x.ampNormClass, data, 'UniformOutput', false), 'Normalized Fluorescence Amplitude', {'immobile', 'confined', 'free', 'directed', 'undetermined'}, 'Normalized intensity (monomer units)');
%plotData(cellfun(@(x) x.ampStatsF20, data, 'UniformOutput', false), '', {'mean', 'first mode mean', 'first mode std', 'first mode fraction', 'number of modes', 'normalized mean'}, '');
%plotData(cellfun(@(x) x.ampStatsL20, data, 'UniformOutput', false), '', {'mean', 'first mode mean', 'first mode std', 'first mode fraction', 'number of modes', 'normalized mean'}, '');
plotData(cellfun(@(x) x.rateMS, data, 'UniformOutput', false), 'Merging and Spliting', {'merging', 'spliting'}, '(per frame)');
%Do only if input specify that this plot be shown. Will cause error if
%data.partitionFrac is not present
if showPartition
    plotData(cellfun(@(x) x.partitionFrac, data, 'UniformOutput', false), 'Partitioning Behavior', {'immobile', 'confined', 'free', 'directed', 'undetermined', 'determined', 'total'}, 'Partition fraction');
end
%% Nested function for plotting and saving
% Splits data structure elements by columns and calls the plotting function
%In other words, converts subData which is cell array of arrays into cell
%array of columns.
%subData            : contains all data to be plotted
%title_Base         : part of plot title that does not change
%title_variable     : part of plot title that does change depending on the
%                     column number
%yLabelName         : used for ylabel of the plot
    function plotData(subData, title_Base, title_Variable, yLabelName)
        nColumns = numel(title_Variable);
        plotData(nColumns) = struct('fitData', [], 'condition', []);
        for iColumns = 1:nColumns
            fitData = plotMultipleSmoothingSpline(outputDirFig, cellfun(@(x) x(:,iColumns), subData, 'UniformOutput', false), times, names, colors, [title_Base '_' title_Variable{iColumns}], yLabelName, smoothingPara);
            plotData(iColumns).fitData = fitData;
            plotData(iColumns).condition = title_Variable{iColumns};
        end
        save([outputDirFig filesep title_Base '.mat'], 'plotData');
    end
%% Save 
end
%% Local Functions





