function [numTrials,tagParam,tagInit] = initCenTagMotion

numTrials = 10;

for i = 1:numTrials
    tagParam(i) = struct('free',1,'chromL',[],'chromS',[],'viscosity',[],...
        'temperature',[],'diffConst',0.1);
end

for i=1:numTrials
    tagInit(i,:) = [0 0 i*0.01]; %in micrometers
end
