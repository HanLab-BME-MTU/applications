function [caseResSummary] = resultsIndTimeCourseMod(ML, saveFile, channels)
%RESULTSINDTIMECOURSE compiles the results of a group of movies making one or multiple timecourse datasets and orders them based on time
%
%SYNOPSIS [caseTimeList,caseResSummary] = resultsIndTimeCourse(ML,caseParam)
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
%
%% DEPRECATION IMMINENT
% mkitti: I am merging this back into resultsIndTimeCourse
%%

%% Input and pre-processing

if nargin < 2
    saveFile = true;
end

if nargin < 3
    channels = [];
end

%get number of movies and number of cases
numMovies = length(ML.movieDataFile_);
numCases = 1;

%reserve memory for individual movie results
resSummary = timeCourseAnalysis.util.emptyResSummaryStruct;

resSummary = repmat(resSummary,numMovies,1);

%define directory for saving results
dir2save = [ML.movieListPath_ filesep 'analysisKJ'];

%% Calculations
nMD = numMovies;
iMD = 0;
%printLength = fprintf(1,'%g/%g MovieData analyzed\n', iMD, nMD);
progressTextMultiple('analyzing MD', nMD)

%go over all movies
for iM = 1 : numMovies
    
    %read in raw results
    %     disp([num2str(iM) '  ' ML.movieDataFile_{iM}]);
    %No need to load MD again if it has already been loaded
    if isempty(ML.movies_{iM})
        MD = MovieData.load(ML.movieDataFile_{iM});
    else
        MD = ML.movies_{iM};
    end
    curChannels = channels(:)';
    if(isempty(curChannels))
        curChannels = 1 : length(MD.channels_);
    else
        curChannels = curChannels(curChannels <= length(MD.channels_));
    end
%     file2savePerMovie = fullfile(dir2save,sprintf([dir2save filesep 'resSummary_movie_%03d.mat'],iM));
    file2savePerMovie = false;
    resSummary(iM,curChannels) = resultsIndTimeCoursePerMovie(MD, file2savePerMovie, curChannels);
end
%fprintf(repmat('\b',1,printLength));

%% Sorting and Saving

%go over each case and put its results together
for iCase = 1 : numCases
    
    %collect and sort results
    caseResSummary = resSummary;
    
    %Saves by default
    if saveFile
    
        %name variables properly for saving
        eval(['timeList_' caseName ' = caseTimeList;'])
        eval(['resSummary_' caseName ' = caseResSummary;'])
        
        %save results
        file2save = fullfile(dir2save,['resSummary_' caseName]); %#ok<NASGU>
        eval(['save(file2save,''timeList_' caseName ''',''resSummary_' caseName ''');']);
        
    end
    
end

