function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK

initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.2,'kTOnTFree',0.22,'addAmpT',0.17,'addWidT',20,'addLenT',0.8,'kTOff',0.0,'kTOnD',0.13,'kDOffFree',8,'addAmpD',4,'addWidD',20,'addLenD',0.6,'kHydrolysis',5.2);

runInfo = struct('maxNumSim',1,'totalTime',600,'simTimeStep',0.01,'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',1);

saveTraj = struct('saveOrNot',0);

saveStats = struct('saveOrNot',0);
