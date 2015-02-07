function [params,dirs] = whInitParamsDirs(params, mainDirname, expname, nFrames)

%% Parameters
if ~isfield(params,'pixelSize') || ~isfield(params,'timePerFrame')
    error('pixelSize and timePerFrame are obligatory parameters');
end

if ~isfield(params,'isDx')
    params.isDx = true;
end


if ~isfield(params,'frameJump')
    params.frameJump = 1;
end

if ~isfield(params,'maxSpeed')
    params.maxSpeed = 90; % um / hr    
end

if ~isfield(params,'fixGlobalMotion')
    params.fixGlobalMotion = false; % correct for bias due to motion
end

if ~isfield(params,'nRois')
    params.nRois = 1; 
end

params.searchRadiusInPixels = ...
    ceil((params.maxSpeed/params.pixelSize)*...
    (params.timePerFrame*params.frameJump/60)); 

params.toMuPerHour = params.pixelSize * 60/(params.timePerFrame*params.frameJump);


if ~isfield(params,'patchSize')
    params.patchSize = ceil(15.0/params.pixelSize); % 15 um in pixels
end

if ~isfield(params,'trajLength')
    params.trajLength = 5;
end

if ~isfield(params,'nBilateralIter')
    params.nBilateralIter = 1;
end

if ~isfield(params,'minClusterSize') % in mu
    params.minClusterSize = 4000;
end

if ~isfield(params,'regionMerginParams')
    params.regionMerginParams.P = 0.03;% small P --> more merging
    params.regionMerginParams.Q = 0.005;% large Q --> more merging (more significant than P)
end

if ~isfield(params,'kymoResolution') % jumps of patchSize
    params.kymoResolution.maxDistMu = 500;
    params.kymoResolution.min = params.patchSize;
    params.kymoResolution.stripSize = params.patchSize;
    params.kymoResolution.max = ceil(params.kymoResolution.maxDistMu/params.pixelSize); % 500 um in pixels    
end

params.strips =  params.kymoResolution.min : params.kymoResolution.stripSize : params.kymoResolution.max;
params.nstrips = length(params.strips);

if ~isfield(params,'maxNFrames')
    params.maxNFrames = 300;%100;
end

if ~isfield(params,'ntime')
    params.nTime = min(nFrames,params.maxNFrames);
end

if ~isfield(params,'always')
    params.always = false;
end


%% Directories

dirs.main = mainDirname;
dirs.dirname = [dirs.main expname];
dirs.expname = expname;

% images
dirs.images = [dirs.dirname '/images/'];

% MF
dirs.mf = [dirs.dirname '/MF/'];
dirs.mfData = [dirs.mf 'mf/'];
dirs.mfScores = [dirs.mf 'scoresVis/'];
dirs.mfBilateral = [dirs.mf 'bilateral/'];

% ROI
dirs.roi = [dirs.dirname '/ROI/'];
dirs.roiData = [dirs.roi 'roi/'];
dirs.roiVis = [dirs.roi 'vis/'];

% Strain rate
dirs.strainRate = [dirs.dirname '/strainRate/'];

% Acceleration
dirs.acceleration = [dirs.dirname '/acceleration/'];

% Coordination
dirs.coordination = [dirs.dirname '/coordination/'];

% kymographs
dirs.kymographs = [dirs.main 'kymographs/'];
dirs.speedKymograph = [dirs.kymographs 'speed/'];
dirs.directionalityKymograph = [dirs.kymographs 'directionality/'];
dirs.strainRateKymograph = [dirs.kymographs 'strainRate/'];
dirs.accelerationKymograph = [dirs.kymographs 'acceleration/'];
dirs.coordinationKymograph = [dirs.kymographs 'coordination/'];

% trajectories
dirs.trajectories = [dirs.main 'trajectories/'];

% Healing rate
dirs.healingRate = [dirs.main 'healingRate/'];
dirs.segmentation = [dirs.main 'segmentation/'];
dirs.plithotaxis = [dirs.main 'plithotaxis/'];

%% Create local directories
if ~exist(dirs.dirname,'dir')
    unix(sprintf('mkdir %s',dirs.dirname));
end

if ~exist(dirs.images,'dir')
    unix(sprintf('mkdir %s',dirs.images));
end

if ~exist(dirs.mf,'dir')
    unix(sprintf('mkdir %s',dirs.mf));
end

if ~exist(dirs.mfData,'dir')
    unix(sprintf('mkdir %s',dirs.mfData));
end

if ~exist(dirs.mfScores,'dir')
    unix(sprintf('mkdir %s',dirs.mfScores));
end

if ~exist(dirs.mfBilateral,'dir')
    unix(sprintf('mkdir %s',dirs.mfBilateral));
end

if ~exist(dirs.roi,'dir')
    unix(sprintf('mkdir %s',dirs.roi));
end

if ~exist(dirs.roiData,'dir')
    unix(sprintf('mkdir %s',dirs.roiData));
end

if ~exist(dirs.roiVis,'dir')
    unix(sprintf('mkdir %s',dirs.roiVis));
end

if ~exist(dirs.strainRate,'dir')
    unix(sprintf('mkdir %s',dirs.strainRate));
end

if ~exist(dirs.acceleration,'dir')
    unix(sprintf('mkdir %s',dirs.acceleration));
end

if ~exist(dirs.coordination,'dir')
    unix(sprintf('mkdir %s',dirs.coordination));
end

%% Global directories
if ~exist(dirs.kymographs,'dir')
    unix(sprintf('mkdir %s',dirs.kymographs));
end

if ~exist(dirs.speedKymograph,'dir')
    unix(sprintf('mkdir %s',dirs.speedKymograph));
end

if ~exist(dirs.directionalityKymograph,'dir')
    unix(sprintf('mkdir %s',dirs.directionalityKymograph));
end

if ~exist(dirs.strainRateKymograph,'dir')
    unix(sprintf('mkdir %s',dirs.strainRateKymograph));
end

if ~exist(dirs.accelerationKymograph,'dir')
    unix(sprintf('mkdir %s',dirs.accelerationKymograph));
end

if ~exist(dirs.coordinationKymograph,'dir')
    unix(sprintf('mkdir %s',dirs.coordinationKymograph));
end

if ~exist(dirs.trajectories,'dir')
    unix(sprintf('mkdir %s',dirs.trajectories));
end

if ~exist(dirs.healingRate,'dir')
    unix(sprintf('mkdir %s',dirs.healingRate));
end

if ~exist(dirs.segmentation,'dir')
    unix(sprintf('mkdir %s',dirs.segmentation));
end

if ~exist(dirs.plithotaxis,'dir')
    unix(sprintf('mkdir %s',dirs.plithotaxis));
end

end