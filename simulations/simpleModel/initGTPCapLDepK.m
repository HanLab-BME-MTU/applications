function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK

initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.5,'kTOnTFree',0.3,...
    'addAmpT',0.2,'addWidT',10,'addLenT',0.85,'kTOff',0.0,'kTOnD',0.2,...
    'kDOffFree',14,'addAmpD',7,'addWidD',10,'addLenD',0.85,'kHydrolysis',10);

saveTraj = struct('saveOrNot',1,'fileName','sampleTraj1');

saveStats = [];
