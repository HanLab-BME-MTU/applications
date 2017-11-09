function [caseResSummary] = resultsIndTimeCourseMod(ML, saveFile, channels, redoPerMovieAnalysis)
%RESULTSINDTIMECOURSEMOD compiles the results of a group of movies making one or multiple timecourse datasets and orders them based on time
%
%SYNOPSIS [caseResSummary] = resultsIndTimeCourseMod(ML,caseParam)
%
%INPUT  ML       : MovieList object containing all movies, either all
%                  belonging to one timecourse.
%        saveFile: Logical that determines if this function saves a file.
%                  The default is 'true'. 0 or 1 instead of true or false
%                  will work.
%       channels : (optional) What channels to use, default:
%                             1:length(MD.channels_)
%    
%OUTPUT 
%       caseResSummary: Structure array storing various results for each
%                       movie, ordered based on their time (corresponding
%                       to caseTimeList). Contains the fields: 
%           .diffSummary         : Diffusion analysis summary, as output
%                                  by summarizeDiffAnRes.
%           .diffCoefMeanPerClass: Mean diffusion coefficient per motion
%                                  class. Order: Immobile, confined, free,
%                                  directed, undetermined.
%           .confRadMeanPerClass : Mean confinement radius per motion
%                                  class. Same order as above.
%           .ampMeanPerClass     : Mean particle amplitude per motion
%                                  class. Rows in same order as above.
%                                  Columns: first for "absolute" amplitude,
%                                  second for amplitude normalized by
%                                  monomer amplitude, as derived from modal
%                                  analysis of particles in last 20 frames
%                                  of each movie.
%           .ampStatsF20         : Amplitude statistics in first 20
%                                  frame of each movie. 
%                                  Order: mean amplitude, first mode
%                                  mean, first mode std, first mode
%                                  fraction, number of modes, mean
%                                  normalized by monomer amplitude (see
%                                  ampMeanPerClass for details).
%           .ampStatsL20         : Same as ampStatsF20, but for last 20
%                                  frames of movie.
%
%Khuloud Jaqaman, March 2015
%modified from resultsIndTimeCourse
%Tae Kim, July 2015
%merged back into resultsIndTimeCourse
%Mark Kittisopikul, January 2016

if nargin < 2
    saveFile = true;
end

if nargin < 3
    channels = [];
end

if nargin < 4
    redoPerMovieAnalysis = true;
end

% mkitti: merged into resultsIndTimeCourse
[~,caseResSummary] = resultsIndTimeCourse(ML,[],saveFile,channels,redoPerMovieAnalysis,'none');
