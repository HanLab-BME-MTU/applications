function [initialState,modelParam,runInfo,saveTraj,saveStats] = initDynInstability

initialState = struct('mtLength0',0.44,'mtState0',1,'free',1,'unitConc',20);

modelParam = struct('minLength',-10000,'kOnElong',0.3,'kOffElong',5,'kOnShrink',0.5,...
    'kOffShrink',11);

runInfo = struct('maxNumSim',1,'totalTime',10000,'simTimeStep',0.01,'timeEps',0.5,...
    'expTimeStep',1,'aveInterval',0.6);

saveTraj = struct('saveOrNot',0);
    
saveStats = struct('saveOrNot',0);
