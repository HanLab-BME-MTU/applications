function [initialState,modelParam,runInfo] = initDynInstability

initialState = struct('mtLength0',0.44,'mtState0',1,'free',1,'unitConc',10);

modelParam = struct('minLength',-10000,'kOnElong',0.6,'kOffElong',4,...
    'kOnShrink',0.2,'kOffShrink',3);

runInfo = struct('maxNumSim',1,'totalTime',100,'simTimeStep',0.001,...
    'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6);

