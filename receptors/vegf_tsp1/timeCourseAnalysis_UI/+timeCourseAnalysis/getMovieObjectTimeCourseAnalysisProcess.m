function [ timeCourseAnalysisProcess ] = getMovieObjectTimeCourseAnalysisProcess( MO , outputArray)
%getMovieObjectTimeCourseAnalysisProcess Obtain TimeCourseAnalysisProcess
%for a MovieObject
%
% See also getMovieDAtaTimeCourseAnalysisProcess

    if(nargin < 2)
        outputArray = true;
    end

    % Actually the MD code is general
timeCourseAnalysisProcess = timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess(MO, outputArray);


end

