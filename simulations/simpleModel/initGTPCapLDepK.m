function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK

initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

runInfo = struct('maxNumSim',1,'totalTime',5000,'simTimeStep',0.004,'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',1);

modelParam = struct('minLength',0.4,'maxLength',0.8,'kTOnTFree',0.96,'addAmpT',0.72,'addWidT',20,'addLenT',0.6,'kTOff',0.0,'kTOnD',0.36,'kDOffFree',120,'addAmpD',108,'addWidD',20,'addLenD',0.6,'kHydrolysis',18);

for i=1:1
    saveTraj(i) = struct('saveOrNot',0,'fileName','tempTraj');
end

saveStats = struct('saveOrNot',0);
