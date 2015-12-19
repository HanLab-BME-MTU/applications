function [p] = UITimeCourseAnalysis(varargin)
%TimeCourseAnalysis of user selected CombinedMovieList objects
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%
%SYNOPSIS [p] = UITimeCourseAnalysis(outputFolder,combinedMovieLists,channelTableData)
%
% INPUT
% outputFolder: absolute path to folder where to save data (optional)
%               as a cell array containing a char
%               default: prompt the user via UI
%
%               Alternatively, the output structure p described below.
%               Other input parameters will be ignored.
% combinedMovieLists: path to combined movie lists as a cell array
%                     (optional)
%                     default: prompt the user via UI
% channelTableData: cell array to insert into the channel table UI
%                   numChannel x 4 cell array
%                   1st column: logical whether channel is selected or not
%                   2nd column: Channel legend
%                   3rd column: Channel name
%                   4th column: Emission wavelength
%                   default: Extract information from first movie list
%
% OUTPUT
% p: parameter structure containing the following fields
%  .outputDir: char. Path to output directory
%  .CML_FullPath: char cell array with paths to combined movie lists
%  .doNewAnalysis: logical. True if time course analysis should rerun
%  .doPartition: logical. True if partitioning analysis should be done
%  .start2zero: logical. True if time should start at zero, false will
%                        align the start time time to zero and there will
%                        be negative values
%  .shiftPlotPositive: logical. True if all timepoints should be positive.
%  .nBootstrp: numermic scalar. Number of permutations to do. default: 100
%  .detectOutliers_k_sigma: numeric scalar. Number of sigma outside of
%                           which that label an outlier
%  .channelTable: cell array containing channel information as described in
%                 channelTableData above
%  .channelNames: Name of channels
%  .channels: Channel numbers to use
%
%Tae H Kim, July 2015
%Mark Kittisopikul, November 2015

%% Initialize
%Progresstext
clear progressTextMultiple;

p = TimeCourseAnalysisConfig(varargin{:});

if(~isempty(p))
    disp(p);
    overallTime = tic;
    timeCourseAnalysis(p.CML_FullPath, p.outputDir, ...
          p ...
          );
%         'doNewAnalysis', p.doNewAnalysis, ...
%         'doPartition', p.doPartition, ...
%         'start2zero', p.start2zero, ...
%         'shiftPlotPositive', p.shiftPlotPositive, ...
%         'channelNames', p.channelTable(:,2), ...
%         'channels',find([p.channelTable{:,1}]), ...
%         'nBootstrp',p.nBootstrp, ...
%         'detectOutliers_k_sigma',p.detectOutliers_k_sigma ...
%         );
    fprintf('Total UITimeCourseAnalysis Elapsed Time %s\n', ...
        char(duration(0,0,toc(overallTime))));
end

