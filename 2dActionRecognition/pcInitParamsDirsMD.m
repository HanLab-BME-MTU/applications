function [params,dirs] = pcInitParamsDirsMD(MD,params)

%% Parameters
if ~isfield(params,'pixelSize') || ~isfield(params,'timePerFrame')
    error('pixelSize and timePerFrame are obligatory parameters');
end

if ~isfield(params,'frameJump')
    params.frameJump = 1;
end

if ~isfield(params,'maxSpeed')
    params.maxSpeed = 60; % um / hr    
end

params.searchRadiusInPixels = ...
    ceil((params.maxSpeed/params.pixelSize)*...
    (params.timePerFrame*params.frameJump/60)); 

params.toMuPerHour = params.pixelSize * 60/(params.timePerFrame*params.frameJump);


if ~isfield(params,'patchSize')
    params.patchSize = ceil(15.0/params.pixelSize); % 15 um in pixels
end

if ~isfield(params,'nTime')
    params.nTime = MD.nFrames_;
end

if ~isfield(params,'always')
    params.always = false;
end

% Fine resolution
if ~isfield(params,'fineResolution')
    params.fineResolution = ceil(3/params.pixelSize); % 3um in pixels (10 pixels)
end

if ~isfield(params,'fineSearchRadius')
    params.fineSearchRadius = ceil(3/params.pixelSize); % 3um in pixels (10 pixels)
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
dirs.roiVis = [dirs.roi 'vis/'];

% Detection of single cells
dirs.detect = [dirs.dirname '/detectCells/'];
dirs.detectData = [dirs.detect '/detections/'];
dirs.detectVis = [dirs.detect '/detectionsVis/'];
dirs.lbp = [dirs.detect '/lbp/'];
dirs.localMorphDynam = [dirs.detect '/localMorphDynam/'];

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

if ~exist(dirs.lbp,'dir')
    unix(sprintf('mkdir %s',dirs.lbp));
end

if ~exist(dirs.localMorphDynam,'dir')
    unix(sprintf('mkdir %s',dirs.localMorphDynam));
end


%% Global Directories
dirs.results = [dirs.dirname '/../results/'];
if ~exist(dirs.results,'dir')
    unix(sprintf('mkdir %s',dirs.results));
end
end