function [params,dirs] = pcProcessTimeLapse(pixelSize, timePerFrame, mainDirname, expname, flags)

params.pixelSize = pixelSize;
params.timePerFrame = timePerFrame;

fname = [mainDirname expname '.tif'];
info = imfinfo(fname);
nFrames = numel(info);

[params,dirs] = pcInitParamsDirs(params, mainDirname, expname,nFrames);

whLocalMotionEstimation(params,dirs);
pcFineMatchingScores(params,dirs);
pcDetectMotionEvents(params,dirs);
% params.always = true;
pcLocalizeMotionEvents(params,dirs);
pcLBP(params,dirs);
pcMorphLocalSpeed(params,dirs);
pcPrepareForTracking(params,dirs);
pcTracking(params,dirs);
end