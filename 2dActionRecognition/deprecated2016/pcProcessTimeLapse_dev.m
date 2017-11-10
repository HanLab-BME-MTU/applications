function [params,dirs] = pcProcessTimeLapse_dev(pixelSize, timePerFrame, mainDirname, expname, flags)

params.pixelSize = pixelSize;
params.timePerFrame = timePerFrame;

fname = [mainDirname expname '.tif'];
info = imfinfo(fname);
nFrames = numel(info);

[params,dirs] = pcInitParamsDirs(params, mainDirname, expname,nFrames);

% whLocalMotionEstimation(params,dirs);
% pcFineMatchingScores(params,dirs);
try
    pcDetectMotionEvents(params,dirs);
catch e
end
% params.always = true;
try
    pcLocalizeMotionEvents(params,dirs);
catch e
end

try    
    pcLBP(params,dirs);
catch e
end

try    
    pcMorphLocalSpeed(params,dirs);
catch e
end

try    
    pcPrepareForTracking(params,dirs);
catch e
end

try
    pcTracking(params,dirs);
catch e
end

try
    pcVisualizeTracking(params,dirs);
catch e
end

end