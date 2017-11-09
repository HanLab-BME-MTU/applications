function [ redoMLanalysis ] = analyzeMDsInParallel( CMLs , doNewAnalysis)
%analyzeMDsInParallel Anlyze all the MovieData objects in a
%CombinedMovieList array

redoMLanalysis = false;
% Analyze all MovieData in parallel first
MLs = [CMLs.movieLists_];
MDs = [MLs.movies_];
% Get TimeCourseAnalysis processes or add them
procs = timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess(MDs);
if(doNewAnalysis)
    % Do per movie part of time course analysis (Alternatively, proc.run)
    MDsummaries = parcellfun_progress(@resultsIndTimeCoursePerMovie,MDs,'UniformOutput',false,'Heading','MovieData analysis');
    % Assign results into processes
    [procs.summary_] = MDsummaries{:};
    cellfun(@save,MDs);
    redoMLanalysis = true;
else
    resSummary = arrayfun(@(proc) proc.summary_,procs,'UniformOutput',false);
    todo = cellfun('isempty',resSummary);
    if(any(todo))
        MDsummaries(todo) = parcellfun_progress(@resultsIndTimeCoursePerMovie,MDs(todo),'UniformOutput',false,'Heading','MovieData analysis');
        % Assign results into processes
        [procs(todo).summary_] = MDsummaries{todo};
        cellfun(@save,MDs);
        redoMLanalysis = true;
    end
end

if(redoMLanalysis)
    % Reset all MovieList level analysis if that needs to be redone
    % Potential efficiency gain: Figure out which ML needs to be reanalyzed
    % by mapping todo back to the ML
    procs = timeCourseAnalysis.getMovieObjectTimeCourseAnalysisProcess(MLs);
    [procs.summary_] = deal([]);
end

end

