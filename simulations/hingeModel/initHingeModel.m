function [hingeParam,hingeInit] = initHingeModel

hingeParam = struct('free',1,'chromL',[],'chromS',[],'viscosity',[],'temperature',[]);

hingeInit = [0 0 0.09]; %in micrometers
