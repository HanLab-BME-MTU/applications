function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK

initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.4,'maxLength',0.8,'kTOnTFree',0.825,'addAmpT',0.6,'addWidT',20,'addLenT',0.6,'kTOff',0.0,'kTOnD',0.3,'kDOffFree',100,'addAmpD',90,'addWidD',20,'addLenD',0.6,'kHydrolysis',15);

runInfo = struct('maxNumSim',5,'totalTime',1000,'simTimeStep',0.005,'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',1);

for i=1:5
    saveTraj(i) = struct('saveOrNot',0,'fileName','tempTraj');
end

saveStats = struct('saveOrNot',0);
