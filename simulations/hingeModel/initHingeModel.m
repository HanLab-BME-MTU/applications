function [hingeParam,hingeInit] = initHingeModel

hingeParam = struct('free',0,'chromL',0.75,'chromS',0.0075,'viscosity',0.2,'temperature',308);

hingeInit = [0 0 0.04];
