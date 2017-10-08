function [] = processTimeLapse(filename,params) % pixelSize, timePerFrame, mainDirname, exp, initParams

[params,dirs] = initParamsDirs(filename,params);

whLocalMotionEstimation(params,dirs);
whTemporalBasedSegmentation(params,dirs);
whCorrectGlobalMotion(params,dirs);
whSegmentationMovie(params,dirs);
whHealingRate(params,dirs); 
whCoordination(params,dirs);
whKymographs(params,dirs);
end