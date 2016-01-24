function [ ] = analyzeMDsInParallel( CMLs , doNewAnalysis)
%analyzeMDsInParallel Anlyze all the MovieData objects in a
%CombinedMovieList array

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
else
    resSummary = arrayfun(@(proc) proc.summary_,procs,'UniformOutput',false);
    todo = cellfun('isempty',resSummary);
    if(any(todo))
        MDsummaries(todo) = parcellfun_progress(@resultsIndTimeCoursePerMovie,MDs(todo),'UniformOutput',false,'Heading','MovieData analysis');
        % Assign results into processes
        [procs(todo).summary_] = MDsummaries{todo};
        cellfun(@save,MDs);
    end
end

end

