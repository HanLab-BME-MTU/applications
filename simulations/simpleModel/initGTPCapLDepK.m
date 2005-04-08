function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK
%INITGTPCAPLDEPK initializes input variables for analyzeMtTrajectory

%If mtGTPCapLDepK alone is called, runInfo and saveStats are not needed.

%See function "mtGTPCapLDepK" for details.
initialState = struct('mtLength0',0.6,'capSize0',3,'unitConc',20);

%See function "analyzeMtTrajectory" for description of variables 
runInfo = struct('maxNumSim',1,'totalTime',2000,'simTimeStep',0.002,'timeEps',0.5,...
    'expTimeStep',0.93,'aveInterval',0.6);

%See function "mtGTPCapLDepK" for meanings and units of variables
modelParam = struct('minLength',0.3,'maxLength',0.9,'kTOnTFree',0.75,'addAmpT',0.72,...
    'addWidT',15,'addLenT',0.62,'kTOff',0.0,'kTOnD',0.4,'kDOffFree',105,'addAmpD',97,...
    'addWidD',15,'addLenD',0.6,'kHydrolysis',14);

%Indicates whether and where generated trajectories trajectories should be saved
saveTraj = struct('saveOrNot',0,'fileName','tempTraj');

%Indicated whether and where trajectory statistics should be saved
saveStats = struct('saveOrNot',0,'fileName','tempStat');
