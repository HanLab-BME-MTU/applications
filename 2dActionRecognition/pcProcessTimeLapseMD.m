function [] = pcProcessTimeLapseMD(MD, timePerFrame, flags)

params.pixelSize = MD.pixelSize_/1000;
params.timePerFrame = timePerFrame;
params.nTime = MD.nFrames_;
params.fixGlobalMotion = false;

[params,dirs] = pcInitParamsDirsMD(MD,params);

whLocalMotionEstimationMD(MD,params,dirs);
pcFineMatchingScoresMD(MD,params,dirs);
pcDetectMotionEventsMD(MD,params,dirs);
pcLocalizeMotionEventsMD(MD,params,dirs);
pcLBPMD(MD,params,dirs);
pcMorphLocalSpeedMD(MD,params,dirs);
end