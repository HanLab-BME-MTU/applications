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

whLocalMotionEstimation(params,dirs);
whTemporalBasedSegmentation(params,dirs);
whCorrectGlobalMotion(params,dirs);
% whSegmentationMovie(params,dirs);
% whHealingRate(params,dirs); % todo: check that this is not affected from frame-frame microscope repeat error

whStrainRate(params,dirs);
whAcceleration(params,dirs);
% whTrajectories(params,dirs);
whCoordination(params,dirs);
whKymographs(params,dirs);

% whSeedsOfPlithotaxis(params,dirs);
% pval = whVectorFlowAccumulation(params,dirs);

% whPlithotaxisTrajectories(params,dirs);
end