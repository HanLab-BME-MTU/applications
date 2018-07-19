function [ analysisTimes, timeLimit, timeLimitIndx ] = getAnalysisTimes( times, timeResolution, inOutFlag )
%getAnalysisTimes Get analysis times
%
% INPUT
% times - from commonInfo
% timeResolution - from params
%
% OUTPUT
% analysisTimes - for commonInfo
% timeLimit - 

times = cellfun(@(x,y) x(y==1), times, inOutFlag,'UniformOutput',false);
timeMax = cellfun(@(x) max(x), times);
timeMin = cellfun(@(x) min(x), times);
%determine the overall range of all conditions (This is useful later when comparing two curves)
timeMaxMax = max(timeMax);
timeMinMin = min(timeMin);
%convert to number divisible by timeResolution
timeMinMin = timeMinMin - mod(timeMinMin, timeResolution);
timeMaxMax = timeMaxMax - mod(timeMaxMax, -timeResolution);
analysisTime_Union = timeMinMin : timeResolution : timeMaxMax;
%Determine the indx and value of timemin and max for each conditions
nConditions = numel(times);
timeLimit = cell(size(times));
timeLimitIndx = cell(size(times));
for iCond = 1:nConditions
    timeLimitIndx{iCond} = [find(analysisTime_Union <= timeMin(iCond), 1, 'last'), find(analysisTime_Union >= timeMax(iCond), 1, 'first')];
    timeLimit{iCond} = [analysisTime_Union(timeLimitIndx{iCond}(1)), analysisTime_Union(timeLimitIndx{iCond}(2))];
end
% timeLimit = reshape(timeLimit,size(data));
%store this in commonInfo
analysisTimes = cellfun(@(x) x(1) : timeResolution : x(2), timeLimit, 'UniformOutput', false);

end

