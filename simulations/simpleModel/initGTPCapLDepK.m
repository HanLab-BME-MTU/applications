function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK
%INITGTPCAPLDEPK initializes input variables for analyzeMtTrajectory

%If mtGTPCapLDepK alone is called, runInfo and saveStats are not needed.

%See function "mtGTPCapLDepK" for details.
initialState.mtLength0 = 0.6;
initialState.capSize0 = 1;
initialState.unitConc = 20;

%See function "mtGTPCapLDepK" for meanings and units of variables
modelParam.kTOnT.lDepend = 0;
modelParam.kTOnT.kmax = 0.05;
modelParam.kTOnT.kmin = 0.2;
modelParam.kTOnT.lmin = 0.2;
modelParam.kTOnT.lmax = 2;

modelParam.kHydr.lDepend = 0;
% modelParam.kHydr.kmax = 1.3;
modelParam.kHydr.kmax = 0.25;
modelParam.kHydr.kmin = 9;
modelParam.kHydr.lmin = 0.2;
modelParam.kHydr.lmax = 2;
modelParam.kHydr.coupled = 1;

modelParam.kTOff.lDepend = 0;
modelParam.kTOff.kmax = 0;
modelParam.kTOff.kmin = [];
modelParam.kTOff.lmin = [];
modelParam.kTOff.lmax = [];

modelParam.kTOnD.lDepend = 0;
modelParam.kTOnD.kmax = 0.02;
modelParam.kTOnD.kmin = 0.35;
modelParam.kTOnD.lmin = 0.2;
modelParam.kTOnD.lmax = 2;

modelParam.kDOff.lDepend = 0;
modelParam.kDOff.kmax = 2;
modelParam.kDOff.kmin = 30;
modelParam.kDOff.lmin = 0.2;
modelParam.kDOff.lmax = 2;

% %See function "analyzeMtTrajectory" for description of variables 
% runInfo = struct('maxNumSim',1,'totalTime',2000,'simTimeStep',0.002,'timeEps',0.5,...
%     'expTimeStep',1,'aveInterval',0.6);
% 
% %Indicates whether and where generated trajectories trajectories should be saved
% saveTraj = struct('saveOrNot',0,'fileName','tempTraj');
% 
% %Indicated whether and where trajectory statistics should be saved
% % saveStats = struct('saveOrNot',0,'fileName','tempStat');
% saveStats = [];

