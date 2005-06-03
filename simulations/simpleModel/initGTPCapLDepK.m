function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK
%INITGTPCAPLDEPK initializes input variables for analyzeMtTrajectory

%If mtGTPCapLDepK alone is called, runInfo and saveStats are not needed.

%See function "mtGTPCapLDepK" for details.
initialState = struct('mtLength0',1,'capSize0',1,'unitConc',20);

%See function "analyzeMtTrajectory" for description of variables 
runInfo = struct('maxNumSim',1,'totalTime',1000,'simTimeStep',0.5,'timeEps',0.9,...
    'expTimeStep',0.01,'aveInterval',1);

%See function "mtGTPCapLDepK" for meanings and units of variables
modelParam = struct('minLength',0.2,'maxLength',2.5,'kTOnTFree',3,'addAmpT',3,...
    'addWidT',100,'addLenT',2.45,'kTOff',0.0,'kTOnD',2.55,'kDOffFree',90,'addAmpD',90,...
    'addWidD',100,'addLenD',0.25,'kHydrolysis',130);

%Indicates whether and where generated trajectories trajectories should be saved
saveTraj = struct('saveOrNot',0,'fileName','tempTraj');

%Indicated whether and where trajectory statistics should be saved
saveStats = struct('saveOrNot',0,'fileName','tempStat');
