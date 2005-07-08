function [initialState,modelParam,runInfo,saveTraj,saveStats] = initDynInstability

initialState = struct('mtLength0',1,'mtState0',1,'free',1,'unitConc',20);

modelParam = struct('minLength',-10000,'kOnElong',0.3,'kOffElong',3,'kOnShrink',0.2,...
    'kOffShrink',7);

% runInfo = struct('maxNumSim',1,'totalTime',1000,'simTimeStep',0.01,'timeEps',0.5,...
%     'expTimeStep',1,'aveInterval',0.6);
runInfo = [];

saveTraj = struct('saveOrNot',0);
    
% saveStats = struct('saveOrNot',0);
saveStats = [];