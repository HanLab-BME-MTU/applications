%deals with individual CML
function [CMLSummary, CMLTime, CMLExtra, startTime] = CMLAnalyze(CML,analysisPara)
    alignEvent = CML.analysisPara_.alignEvent;
    %[CMLSummary, CMLTime, CMLExtra] = arrayfun(@(x) MLAnalyze(x, alignEvent), CML.movieLists_, 'UniformOutput', false);
    nML = numel(CML.movieLists_);
    for iML = nML:-1:1
        [CMLSummary{iML}, CMLTime{iML}, CMLExtra{iML}, startTime(iML)] = timeCourseAnalysis.MLAnalyze(CML.movieLists_(iML), alignEvent,analysisPara);
    end
    CMLSummary = vertcat(CMLSummary{:});
    CMLTime = vertcat(CMLTime{:});
    CMLExtra = vertcat(CMLExtra{:});
    startTime = mean(startTime);
end
