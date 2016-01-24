%% Time Course Analysis (ML-level)
function [MLSummary, MLTime, MLExtra, startTime] = MLAnalyze(ML, alignEvent,analysisPara)
%Basic Analysis-------------------------------
%new analysis if doNewAnalysis is true or analysis has not been
%done yet
% ML.sanityCheck;
TCAPIndx = ML.getProcessIndex('TimeCourseAnalysisProcess');
%progressText
% progressTextMultiple('part 1', 2);
%(I'm not exactly sure what resultsIndTimeCourseMod does)
%It is used like blackbox that does the basic analysis
if isempty(TCAPIndx) || isempty(ML.processes_{TCAPIndx}.summary_)
    MLSummary = resultsIndTimeCourseMod(ML, false, analysisPara.channels , analysisPara.doNewAnalysis);
    ML.addProcess(TimeCourseAnalysisProcess(ML));
    TCAPIndx = ML.getProcessIndex('TimeCourseAnalysisProcess');
    ML.processes_{TCAPIndx}.setSummary(MLSummary);
    ML.save;
elseif analysisPara.doNewAnalysis
    MLSummary = resultsIndTimeCourseMod(ML, false,analysisPara.channels, true);
    ML.processes_{TCAPIndx}.setSummary(MLSummary);
    ML.save;
else
    MLSummary = ML.processes_{TCAPIndx}.summary_;
end
% progressTextMultiple('part 2');
%Time Analysis---------------------------------
timeProcIndx = ML.getProcessIndex('TimePoints');
alignIndx = ML.processes_{timeProcIndx}.getIndex(alignEvent);
offSet = datenum(ML.processes_{timeProcIndx}.times_{alignIndx});
MLTime = cellfun(@(x) (datenum(x.acquisitionDate_) - offSet) .* 1440, ML.movies_, 'UniformOutput', false);
MLTime = vertcat(MLTime{:});
%get relative startTime
startIndx = ML.processes_{timeProcIndx}.getIndex('start');
startTime = (datenum(ML.processes_{timeProcIndx}.times_{startIndx}) - offSet) * 1440;
%Conditional (Extra) Analysis-------------------------
MLExtra = cellfun(@timeCourseAnalysis.MDAnalyze, ML.movies_, 'UniformOutput', false);
MLExtra = vertcat(MLExtra{:});
%--------------------TrackPartitioning analysis---------------
if analysisPara.doPartition
    %Column 1 : immobile
    %Column 2 : confined
    %Column 3 : free
    %Column 4 : directional
    %Column 5 : undetermined
    FN = {'immobile', 'confined', 'free', 'directed', 'undetermined'};
    nFN = numel(FN);
    [result_pCoef, ~] = partitionCoef(ML);
    nMD = numel(ML.movies_);
    for iMD = 1:nMD
        for iFN = 1:nFN
            MLExtra(iMD,1).chemEnergy(iFN) = log(result_pCoef.partCoef(iMD).(FN{iFN}));
            MLExtra(iMD,1).locFreq(iFN) = result_pCoef.locFreq(iMD).(FN{iFN});
            MLExtra(iMD,1).delocFreq(iFN) = result_pCoef.delocFreq(iMD).(FN{iFN});
            MLExtra(iMD,1).eqCond(iFN) = result_pCoef.eqCond(iMD).(FN{iFN});
        end
    end
end
% progressTextMultiple();
%---------------------
%progress text
progressTextMultiple();
end