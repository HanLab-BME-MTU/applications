function [numTrials,hingeParam,hingeInit] = initHingeModel

numTrials = 10;

for i = 1:numTrials
    hingeParam(i) = struct('free',1,'chromL',[],'chromS',[],'viscosity',[],...
        'temperature',[],'diffConst',0.1);
end

for i=1:numTrials
    hingeInit(i,:) = [0 0 i*0.01]; %in micrometers
end
