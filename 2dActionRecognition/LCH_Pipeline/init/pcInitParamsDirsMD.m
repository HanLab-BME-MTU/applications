function [params,dirs] = pcInitParamsDirsMD(MD,curFname,curTask,params)

%% Parameters
if ~isfield(params,'pixelSize') || ~isfield(params,'timePerFrame')
    error('pixelSize and timePerFrame are obligatory parameters');
end

if ~isfield(params,'curFname')
    params.curFname = curFname;
end

if ~isfield(params,'curTask')
    params.curTask = curTask;
end

if ~isfield(params,'frameJump')
    params.frameJump = 1;
end

if ~isfield(params,'maxSpeed')
    params.maxSpeed = 90; % um / hr    
end

params.searchRadiusInPixels = ...
    ceil((params.maxSpeed/params.pixelSize)*...
    (params.timePerFrame*params.frameJump/60)); 

params.toMuPerHour = params.pixelSize * 60/(params.timePerFrame*params.frameJump);


if ~isfield(params,'patchSize')
    params.patchSize = ceil(10.0/params.pixelSize); % 10 um in pixels
end

if ~isfield(params,'nTime')
    params.nTime = MD.nFrames_;
end

if ~isfield(params,'sTime')
    params.sTime = 1;
end

if ~isfield(params,'always')
    params.always = false;
end

if ~isfield(params,'deepDebug')
    params.always = false;
end

% Fine resolution
if ~isfield(params,'fineResolution')
    params.fineResolution = params.searchRadiusInPixels;%ceil(3/params.pixelSize); % 3um in pixels (10 pixels)
end

if ~isfield(params,'fineSearchRadius')
    params.fineSearchRadius = params.searchRadiusInPixels;%ceil(3/params.pixelSize); % 3um in pixels (10 pixels)
end

% Single cell parameters
if ~isfield(params,'cellMinArea') % in pixels
    params.cellMinArea = (5/params.pixelSize)^2;
end
if ~isfield(params,'cellMaxArea') % in pixels
    params.cellMaxArea = (35/params.pixelSize)^2;
end
if ~isfield(params,'singleCellBB4Vis')
    params.singleCellBB4Vis = round(100 ./ params.pixelSize);
end

% what percentile to use as "background" in cross correlation matching
% score 
if ~isfield(params,'detectionPPVicinityPrctile')
    params.detectionPPVicinityPrctile = 10; % percentile
end

% how many farmes to look back and forth to filter out a detection because
% the same hit was filtered earlier
if ~isfield(params,'detectionPPFilterFPTime')
    params.detectionPPFilterFPTime = 5; % # frames
end

% vicinity (in pixels) around a detection to calculate statistics at
if ~isfield(params,'detectionPPVicinityPixels')
    %     params.detectionPPVicinityPixels = 30 * params.pixelSize; % 30um in pixels
    params.detectionPPVicinityPatch = ceil(80 ./ params.pixelSize / params.patchSize); % 80um in patches
end

if ~isfield(params,'diffDetectionBackgroundScoreTH')
    params.diffDetectionBackgroundScoreTH = 0.1; % based on the data from Analysis/metaAnalysis
end

if ~isfield(params, 'minLengthTH')
    params.minLengthTH = 120; %45
end

if ~isfield(params, 'minLengthTHLong')
    params.minLengthTHLong = 240;
end

% Radius for a single cell for the LBP
if ~isfield(params, 'FOVRadius')
    params.FOVRadius = round(35/params.pixelSize); % 35 um in pixels
end

% Radius for Otsu cell segmentaion
if ~isfield(params, 'RoiRadius')
    params.RoiRadius = round(120/params.pixelSize); % 120 um in pixels
end

if ~isfield(params, 'nScales')
    params.nScales = 4; % 4 scales (1 to 1/8)
end

if ~isfield(params, 'scales')
    params.scales = 1.0./2.^((1:params.nScales)-1);
end

%% Local Directories
dirs.dirname = MD.outputDirectory_;

% images
dirs.images = [dirs.dirname '/images/'];

% MF
dirs.mf = [dirs.dirname '/MF/'];
dirs.mfData = [dirs.mf 'mf/'];
dirs.mfScores = [dirs.mf 'scoresVis/'];

% Fine resolution
dirs.fineRes = [dirs.mf 'fineResolution/'];
dirs.fineResScores = [dirs.fineRes 'scores/'];
dirs.fineResScoresVis = [dirs.fineRes 'scoresVis/'];

% ROI
dirs.roi = [dirs.dirname '/ROI/'];
dirs.roiData = [dirs.roi 'roi/'];
dirs.roiLever = [dirs.roi 'lever/'];
dirs.roiVis = [dirs.roi 'vis/'];

% Detection of single cells
dirs.detect = [dirs.dirname '/detectCells/'];
dirs.detectData = [dirs.detect '/detections/'];
dirs.detectVis = [dirs.detect '/detectionsVis/'];

dirs.detectPPData = [dirs.detect '/detectionsPP/'];
dirs.detectPPVis = [dirs.detect '/detectionsPPVis/'];

dirs.tracking = [dirs.dirname '/tracking/'];

dirs.lbp = [dirs.dirname '/lbp/'];

dirs.lbpDt = [dirs.dirname '/lbpDt/'];

% if exist([dirs.detect '/lbp/'],'dir')
%    unix(sprintf('rm -rf %s',[dirs.detect '/lbp/']));
% end
% 
% if exist([dirs.detect '/localMorphDynam/'],'dir')
%     unix(sprintf('rm -rf %s',[dirs.detect '/localMorphDynam/']));
% end



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

if ~exist(dirs.fineRes,'dir')
    unix(sprintf('mkdir %s',dirs.fineRes));
end

if ~exist(dirs.fineResScores,'dir')
    unix(sprintf('mkdir %s',dirs.fineResScores));
end

if ~exist(dirs.fineResScoresVis,'dir')
    unix(sprintf('mkdir %s',dirs.fineResScoresVis));
end

if ~exist(dirs.roi,'dir')
    unix(sprintf('mkdir %s',dirs.roi));
end

if ~exist(dirs.roiData,'dir')
    unix(sprintf('mkdir %s',dirs.roiData));
end

if ~exist(dirs.roiLever,'dir')
    unix(sprintf('mkdir %s',dirs.roiLever));
end

if ~exist(dirs.roiVis,'dir')
    unix(sprintf('mkdir %s',dirs.roiVis));
end

if ~exist(dirs.detect,'dir')
    unix(sprintf('mkdir %s',dirs.detect));
end

if ~exist(dirs.detectData,'dir')
    unix(sprintf('mkdir %s',dirs.detectData));
end

if ~exist(dirs.detectVis,'dir')
    unix(sprintf('mkdir %s',dirs.detectVis));
end

if ~exist(dirs.detectPPData,'dir')
    unix(sprintf('mkdir %s',dirs.detectPPData));
end

if ~exist(dirs.detectPPVis,'dir')
    unix(sprintf('mkdir %s',dirs.detectPPVis));
end

if ~exist(dirs.tracking,'dir')
    unix(sprintf('mkdir %s',dirs.tracking));
end

if ~exist(dirs.lbp,'dir')
    unix(sprintf('mkdir %s',dirs.lbp));
end

if ~exist(dirs.lbpDt,'dir')
    unix(sprintf('mkdir %s',dirs.lbpDt));
end

% if ~exist(dirs.localMorphDynam,'dir')
%     unix(sprintf('mkdir %s',dirs.localMorphDynam));
% end


%% Global folders
tmp = strsplit(dirs.dirname,filesep);
dirs.expname = tmp{end};

%% Global Directories
dirs.trackingVis = [dirs.dirname '/../../../Movies/trackingMovies/'];
if ~exist(dirs.trackingVis,'dir')
    unix(sprintf('mkdir %s',dirs.trackingVis));
end

dirs.lever = [dirs.dirname '/../../../LEVER/masks/'];

end