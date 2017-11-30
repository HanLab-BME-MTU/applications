function [] = pcProcessTimeLapseMD(MD, curFname, curTask, timePerFrame, flags)

params.pixelSize = MD.pixelSize_/1000;

% Hacking
if MD.pixelSize_< 325
    params.pixelSize = 325/1000;
end

params.timePerFrame = timePerFrame;
params.nTime = min(310,MD.nFrames_);%MD.nFrames_;
params.sTime = 60;%MD.nFrames_;
params.fixGlobalMotion = false;
params.deepDebug = true;

[params,dirs] = pcInitParamsDirsMD(MD,curFname,curTask,params);

% params.always = true;

% pcMD2tifs(MD,params,dirs);

whLocalMotionEstimationMD(MD,params,dirs);

pcStageLocationError(MD,params,dirs);

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
% pcTrackingMovieLong(MD, params, dirs); % Also single cell definition - long trajectories

% params.always = tmpAlways;


% try
% pcVisualizeTracking(MD,params,dirs);
% catch ee
%     return;
% end

%% Single cell segmentaion
params.always = true;
pcCellRoiMD(MD,params,dirs);
pcCellRoiLeverMD(MD,params,dirs);
params.always = false;

%% Feature extraction - per cell!
% params.always = true;
tic;
pcLBPMD(MD,params,dirs);% radius = 35um
pcLBPdtMD(MD,params,dirs);% radius = 35um
% pcLBPMD_LEVER(MD,params,dirs);% radius = 35um
toc

pcSingleCellMovies(MD,params,dirs);
pcSingleCell_dLBP(params,dirs);
pcSingleCell_dLBP_Plasticity(params,dirs)

end