function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK

initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.0,'kTOnTFree',1.2,'addAmpT',0.9,'addWidT',100,'addLenT',0.5,'kTOff',0.0,'kTOnD',0.4,'kDOffFree',40,'addAmpD',13,'addWidD',100,'addLenD',0.4,'kHydrolysis',15);

runInfo = struct('maxNumSim',1,'totalTime',1000,'simTimeStep',0.01,'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',1);

saveTraj = struct('saveOrNot',1,'fileName','tempTraj');

saveStats = struct('saveOrNot',0);
