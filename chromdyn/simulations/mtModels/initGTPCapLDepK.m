function [initialState,modelParam,runInfo,saveTraj,saveStats] = initGTPCapLDepK
%INITGTPCAPLDEPK initializes input variables for analyzeMtTrajectory

%If mtGTPCapLDepK alone is called, runInfo and saveStats are not needed.

%See function "mtGTPCapLDepK" for details.
initialState.mtLength0 = 3;
initialState.capSize0 = 1;
initialState.unitConc = 10;

%See function "mtGTPCapLDepK" for meanings and units of variables
modelParam.kTOnT.lDepend = -1;
modelParam.kTOnT.kmax = 0.05;
modelParam.kTOnT.kmin = 0.04;
modelParam.kTOnT.lmin = 2.5;
modelParam.kTOnT.lmax = 3.5;

modelParam.kHydr.lDepend = 0;
modelParam.kHydr.kmax = 0.25;
modelParam.kHydr.kmin = [];
modelParam.kHydr.lmin = [];
modelParam.kHydr.lmax = [];
modelParam.kHydr.coupled = 1;

modelParam.kTOff.lDepend = 0;
modelParam.kTOff.kmax = 0;
modelParam.kTOff.kmin = [];
modelParam.kTOff.lmin = [];
modelParam.kTOff.lmax = [];

modelParam.kTOnD.lDepend = -1;
modelParam.kTOnD.kmax = 0.02;
modelParam.kTOnD.kmin = 0.015;
modelParam.kTOnD.lmin = 2.5;
modelParam.kTOnD.lmax = 3.5;

modelParam.kDOff.lDepend = 0;
modelParam.kDOff.kmax = 1.5;
modelParam.kDOff.kmin = [];
modelParam.kDOff.lmin = [];
modelParam.kDOff.lmax = [];

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

