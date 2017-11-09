function [pval] = whProcessTimeLapse(pixelSize, timePerFrame, mainDirname, exp, initParams)

if nargin == 5
    params = initParams;
end

params.pixelSize = pixelSize;
params.timePerFrame = timePerFrame;


% fname = [mainDirname exp.name exp.ext];
% nFrames = exp.nFrames;

if isfield(exp,'name')
    [params,dirs] = whInitParamsDirs(params, mainDirname, exp.name, exp.nFrames);
else
    nFrames = 100;
    [params,dirs] = whInitParamsDirs(params, mainDirname, exp, nFrames);
end

% Just for debug!
% params.nTime = 30;

whLocalMotionEstimation(params,dirs);
 
whTemporalBasedSegmentation(params,dirs);
whCorrectGlobalMotion(params,dirs);
% % % whVectorFieldsVisualization(params,dirs);
% % params.always = false;
whSegmentationMovie(params,dirs);
whHealingRate(params,dirs); % todo: check that this is not affected from frame-frame microscope repeat error
% 
% % whStrainRate(params,dirs);
% % whAcceleration(params,dirs);
% % % whTrajectories(params,dirs);
% % params.always = true;
whCoordination(params,dirs);
whKymographs(params,dirs);

whKymographsStd(params,dirs);


whSeedsOfPlithotaxis(params,dirs);
% pval = whVectorFlowAccumulation(params,dirs); % permission denied??

% whPlithotaxisTrajectories(params,dirs);
end