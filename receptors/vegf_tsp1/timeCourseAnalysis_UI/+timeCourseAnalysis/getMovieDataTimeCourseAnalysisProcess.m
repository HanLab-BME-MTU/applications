function [ timeCourseAnalysisProcess ] = getMovieDataTimeCourseAnalysisProcess( MD , outputArray)
%getMovieDataTimeCourseAnalysisProcess Obtains a TimeCourseAnalysisProcess
%associated with a MovieData
% MD : MovieData , or an array of MovieData, or a cell array
% outputArray : true if an array of processes should be output (default)
%               or false if a cell array match the size of MD should be
%               output
    if(nargin < 2)
        outputArray = true;
    end

    % Deal with cell or arrays of MovieData
    if(iscell(MD))
        timeCourseAnalysisProcess = cellfun( ...
            @timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess, ...
            MD, ...
            'UniformOutput',false);
        if(outputArray)
            timeCourseAnalysisProcess = [timeCourseAnalysisProcess{:}];
        end
        return;
    elseif(~isscalar(MD))
        timeCourseAnalysisProcess = arrayfun( ...
            @timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess, ...
            MD, ...
            'UniformOutput',false);
        if(outputArray)
            timeCourseAnalysisProcess = [timeCourseAnalysisProcess{:}];
        end
        return;
    end

    % Check if TimeCourseAnalysisProcess already exists
    idx = MD.getProcessIndex('TimeCourseAnalysisProcess');
    if(isempty(idx))
        % If not, create a new one
        timeCourseAnalysisProcess = TimeCourseAnalysisProcess(MD);
        MD.addProcess(timeCourseAnalysisProcess);
    else
        timeCourseAnalysisProcess = MD.getProcess(idx);
    end
    timeCourseAnalysisProcess.sanityCheck();
end

