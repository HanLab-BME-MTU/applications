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

% params.always = true;
pcSetSingleCellTrajectories(params,dirs);
% params.always = false;

%% movies
% pcDetectionMovie(MD, params, dirs);

% params.always = true;

% is pcTracking was performed than run again pcTrackingMovie (for single cell labeling)
% tmpAlways = params.always;
% params.always = doneTracking || tmpAlways;
% params.always = true;

pcTrackingMovie(MD, params, dirs); % (Once did the single cell definition)

% pcTrackingMovieLong(MD, params, dirs); % Also single cell definition - long trajectories

% params.always = tmpAlways;


% try
% pcVisualizeTracking(MD,params,dirs);
% catch ee
%     return;
% end

%% Single cell segmentaion
pcCellRoiMD(MD,params,dirs);
pcCellRoiLeverMD(MD,params,dirs);


%% Feature extraction - per cell!
% params.always = true;
% tic;
pcLBPMD(MD,params,dirs);% radius = 35um
pcLBPdtMD(MD,params,dirs);% radius = 35um

params.always = true;
curDir = pwd;
cd '/home2/azaritsky/code/extern/hctsa';

% patch to check if the problem is in a specific node
for iInstall = 1 : 100
    try
        install;
    catch ee
        warning(ee.message);
        continue;
    end
    break;
end
if iInstall == 100
    install;
end

pcLBPMD_LEVER(MD,params,dirs);% radius = 35um
pcShapeMD_LEVER(MD,params,dirs);% radius = 35um
eval(['cd ' curDir]);

% toc
params.always = true;

% pcSingleCellMovies(MD,params,dirs);
% pcSingleCell_dLBP(params,dirs);
% pcSingleCell_dLBP_Plasticity(params,dirs);
pcSingleCell_LEVER_LBP_SHAPE(params,dirs);

end