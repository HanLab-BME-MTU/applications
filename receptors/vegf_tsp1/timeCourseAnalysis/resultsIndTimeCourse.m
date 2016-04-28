function [caseTimeList,caseResSummary] = resultsIndTimeCourse(ML,caseParam,saveFile, channels, redoPerMovieAnalysis,parallel)
%RESULTSINDTIMECOURSE compiles the results of a group of movies making one or multiple timecourse datasets and orders them based on time
%
%SYNOPSIS [caseTimeList,caseResSummary] = resultsIndTimeCourse(ML,caseParam)
%
%INPUT  ML       : MovieList object containing all movies, either all
%                  belonging to one timecourse, or potentially belonging
%                  to multiple timecourses. The number of timecourses is
%                  indicated in the next variable.
%       caseParam: Structure array with number of entries = number of
%                  timecourses among which the movies will be divided. For
%                  each time course, it contains the following fields:
%           .indx    : Indices of movies in ML belonging to this
%                      timecourse.
%           .timeList: The time point of each movie, with order following
%                      the movie order in ML.
%           .name    : Case name, to be used in naming output variables.
%           .indx0min: Index of movie to be considered at time 0, e.g. when
%                      a ligand is added. If > 1, movies before will have a
%                      negative relative time, movies after wil have a
%                      positive relative time. If 1, relative time and
%                      absolute time are the same.
%        saveFile: Logical that determines if this function saves a file.
%                  The default is 'true'. 0 or 1 instead of true or false
%                  will work.
%        channels: Index of which channels to analyze (numeric array)
%                  Empty (default) will result in each channel being analyzed
%        redoAnalysis: Logical that determines if analysis should be redone
%                  per movie if prior analysis per movie exists.
%                  The default is true. 0 or 1 instead of true or false
%                  will work.
%        parallel    : char that is either
%                      'none' (default): Run in a single thread
%                      'cluster'       : Run using parallel pool
%                      'slurm'         : Queue on cluster using sbatch
%                      'load_only'     : Load saved data only
%    
%OUTPUT caseTimeList  : 2-column vector indicating the movie times. Column 1
%                       shows absolute time, Column 2 shows relative time.
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

%% Input and pre-processing

if nargin < 1
    error('resultsIndTimeCourse: Too few input arguments')
end
% For a single parameter, see end of file
if nargin < 2
    caseParam = [];
end
if nargin < 3
    saveFile = true;
end
if nargin < 4
    channels = [];
end
if nargin < 5
    redoPerMovieAnalysis = true;
end
if nargin < 6
    parallel = 'none';
end

%get number of movies
% numMovies = length(ML.movieDataFile_);



%do not need to reserve memory for individual movie results
% since cellfun takes care of this


%define directory for saving results
dir2save = [ML.movieListPath_ filesep 'analysisKJ'];
if(~exist(dir2save,'dir'))
    mkdir(dir2save);
end

%% Calculations
% Most of the analysis is done before getting here by running the MD-level
% processes. All we need to do is retrieve that 
procs = timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess(ML.movies_,false);
resSummary = cellfun(@(proc) proc.summary_,procs,'UniformOutput',false);
outFilePaths = cellfun(@(proc) proc.outFilePaths_,procs,'UniformOutput',false);
file2savePerMovie = strcat(outFilePaths,[filesep 'resSummary.mat']);
if(redoPerMovieAnalysis)
    todo = true(size(ML.movieDataFile_));
else
    todo = cellfun('isempty',resSummary);
end

switch(parallel)
    % mkitti: Moved most of the body to resultsIndTimeCoursePerMovie
    case 'none'
        progressTextMultiple('analyzing MD', length(ML.movies_(todo)));
        resSummary(todo) = cellfun(@(dataFile,file2save) resultsIndTimeCoursePerMovie(dataFile,file2save,channels), ...
            ML.movies_(todo), ...
            file2savePerMovie(todo), ...
            'UniformOutput',false);
%         cellfun(@progressTextMultiple,resSummary(todo));
    case 'cluster'
        movieFiles = distributed(ML.movieDataFile_(todo));
        file2savePerMovie = distributed(file2savePerMovie(todo));
        resSummary(todo) = gather(cellfun(@(dataFile,file2save) resultsIndTimeCoursePerMovie(dataFile,file2save,channels), ...
            movieFiles, ...
            file2savePerMovie, ...
            'UniformOutput',false));
    case 'batch'
        cluster = parcluster('nucleus2015a');
        jobs = cellfun(...
            @(dataFile,file2save) ...
                batch(cluster,@resultsIndTimeCoursePerMovie,1, ...
                    {dataFile,file2save,channels}, ...
                    'Pool',31, ...
                    'AutoAttachFiles',false), ...
            ML.movieDataFile_(todo), ...
            file2savePerMovie(todo));
        cellfun(@submit,jobs);
        cellfun(@wait,jobs);
        resSummary(todo) = cellfun(@fetchAndUnwrap,jobs,'UniformOutput',false);
    case 'createJob'
        cluster = parcluster('nucleus2015a');
        job = createJob(cluster,'AutoAttachFiles',false);
        tasks = cellfun(...
            @(dataFile,file2save) ...
                createTask(@resultsIndTimeCoursePerMovie,1, ...
                    {dataFile,file2save,channels}), ...
            ML.movieDataFile_(todo), ...
            file2savePerMovie(todo));
        if(~isempty(tasks))
                submit(job);
        end
        wait(job);
        resSummary(todo) = fetchOutputs(job);
    case 'parcellfun_progress'
        resSummary(todo) = parcellfun_progress(@(dataFile,file2save) resultsIndTimeCoursePerMovie(dataFile,file2save,channels), ...
            ML.movieDataFile_(todo), ...
            file2savePerMovie(todo), ...
            'UniformOutput',false);
    case 'load_only'
        resSummary(todo) = cellfun(@pauseUntilExistThenLoad,ML.movieDataFile_(todo),'UniformOutput',false);
    otherwise
        error(['Invalid parallel parameter: ' parallel]);
end

% Convert resSummary into a numMovies x max(channels) struct array
resSummary = vertcat(resSummary{:});

%% Sorting and Saving

%go over each case and put its results together
[caseTimeList,caseResSummary] = resultsIndTimeCourseSortAndSave(resSummary,caseParam,saveFile);

end
function resSummary = loadResSummary(file)
    s = load(file);
    resSummary = s.resSummary;
end
function resSummary = pauseUntilExistThenLoad(file)
    % Block until file exists
    while(~exist(file,'file'))
        % Pause for 5 seconds if analysis file does not exist
        disp(['Waiting 5 seconds for movie ' num2str(iM) ' out of ' num2str(numMovies)]);
        pause(5);
    end
    resSummary = loadResSummary(file);
end
function resSummary = fetchAndUnwrap(job)
    c = fetchOutputs(job);
    resSummary = c{1};
end