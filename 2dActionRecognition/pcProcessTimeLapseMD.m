function [] = pcProcessTimeLapseMD(MD, curFname, curTask, timePerFrame, flags)

params.pixelSize = MD.pixelSize_/1000;

% Hacking
if MD.pixelSize_< 325
    params.pixelSize = 325/1000;
end

params.timePerFrame = timePerFrame;
params.nTime = MD.nFrames_;
params.fixGlobalMotion = false;
params.deepDebug = true;

[params,dirs] = pcInitParamsDirsMD(MD,curFname,curTask,params);

% params.always = true;

whLocalMotionEstimationMD(MD,params,dirs);

pcDetectionMD(MD,params,dirs);

% params.always = true;

pcDetectionPPCalcLocalMatchScoreMD(params,dirs); % todo:retest the filtered (magenta) detections?

pcDetectionPPFilterLocalMatchScore(MD,params,dirs);

pcDetectionPPBackTrack(MD,params,dirs);

% params.always = true;
doneTracking = pcTracking(params,dirs);

%% movies
% pcDetectionMovie(MD, params, dirs);

% params.always = true;

% is pcTracking was performed than run again pcTrackingMovie (for single cell labeling)
% tmpAlways = params.always;
% params.always = doneTracking || tmpAlways;
% params.always = true;
pcTrackingMovie(MD, params, dirs); % This function also does the single cell definition!
pcTrackingMovieLong(MD, params, dirs); % Also single cell definition - long trajectories

% params.always = tmpAlways;


% try
% pcVisualizeTracking(MD,params,dirs);
% catch ee
%     return;
% end

% params.always = true;
pcCellRoiMD(MD,params,dirs);

params.always = true;
tic;
pcLBPMD(MD,params,dirs);% radius = 35um
toc

pcSingleCellMovies(MD,params,dirs);
pcSingleCell_dLBP(params,dirs);

% TO IMPLEMENT NEXT:
% pcSingleCell_MIP(params,dirs);
return;

whCorrectGlobalMotion(params,dirs); % should be here? do we have infrastructure for that?
pcMorphLocalSpeedMD(MD,params,dirs);
end