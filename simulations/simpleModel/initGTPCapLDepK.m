function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK

initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.5,'kTOnTFree',1.2,'addAmpT',0.4,'addWidT',10,'addLenT',0.85,'kTOff',0.0,'kTOnD',0.4,'kDOffFree',45,'addAmpD',21,'addWidD',10,'addLenD',0.85,'kHydrolysis',25);

runInfo = struct('maxNumSim',1,'totalTime',1000,'simTimeStep',0.005,'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',1);

saveTraj = struct('saveOrNot',0);

saveStats = struct('saveOrNot',0);
