function [hingeParam,hingeInit] = initHingeModel

hingeParam = struct('free',0,'chromL',1,'chromS',0.2,'viscosity',0.2,'temperature',308);

%hingeParam = struct('free',1);

hingeInit = [0 0 0.05]; %in micrometers
