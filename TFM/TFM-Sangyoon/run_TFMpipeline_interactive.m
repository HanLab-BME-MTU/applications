%% run_TFMpipeline_safe2.m
%
% SAFE version of run_TFMpipeline_interactive.m
%
% Pipeline:
%   1) Interactive param prompt (CLI)
%   2) Pick movieList.mat files using old-style uigetfile, not uifigure
%   3) Segmentation: ThresholdProcess on CELL channel
%   4) MaskRefinementProcess with default settings
%   5) TFMPackage setup: bead channel for drift correction + displacement
%   6) Run tfmRun movie-by-movie with progress
%
% Channel conventions:
%   iBeadChan  = bead channel -> EfficientSubpixelRegistration + displacement
%   segChan    = cell channel -> SegmentationPackage / ThresholdProcess
%
% Sangyoon Han

clear; clc;

%% ===== DEFAULTS =====
useParallel = false; %#ok<NASGU> % tfmRun is called movie-by-movie below

defaults = struct();
defaults.iBeadChan    = 2;      % bead channel (for TFMPackage)
defaults.segChan      = 1;      % cell channel (for Segmentation)
defaults.YoungModulus = 2700;   % Pa
defaults.PoissonRatio = 0.49;
defaults.regParam     = 0.005;
defaults.segModeIdx   = 2;      % 1=MinMax, 2=Otsu, 3=Rosin, 4=Gradient
defaults.segGaussSigma= 2;      % pixels (0 = no blur)
defaults.forceRerunSeg= false;
defaults.minCorLength = 17;     % template window half-width (pixels)
defaults.maxFlowSpeed = 20;     % max displacement (pixels per frame)
defaults.addNonLocMaxBeads = false;  % add non-local-max beads (default off)

%% ===== INTERACTIVE PROMPT =====
tfmParams = promptTFMSegParams(defaults);

%% ===== PICK MOVIELISTS =====
% SAFE MODE: avoid uifigure/uilistbox hang.
% Option A: paste one or more movieList.mat paths below and skip GUI.
% Example:
% MLpaths = {'/mnt/nas/path/to/movieList.mat'};
MLpaths = {};

% Option B: if MLpaths is empty, use old-style uigetfile.
% This avoids the uifigure listbox that caused MATLAB/desktop freezing.
if isempty(MLpaths)
    MLpaths = pickMovieListsSafeUigetfile();
end
assert(~isempty(MLpaths), 'No movieList.mat selected.');

%% ===== Ensure SegmentationPackage on path =====
ensureSegmentationPackageOnPath();

%% ===== MAIN LOOP =====
nML = numel(MLpaths);
for ci = 1:nML

    mlPath = char(MLpaths{ci});
    tML = tic;
    fprintf('\n====================================================\n');
    fprintf('[%s] MovieList %d/%d\n  %s\n', datestr(now,'HH:MM:SS'), ci, nML, mlPath);
    fprintf('  beadChan=%d | segChan=%d | E=%.3g Pa | reg=%.3g\n', ...
        tfmParams.iBeadChan, tfmParams.segChan, ...
        tfmParams.YoungModulus, tfmParams.regParam);
    fprintf('====================================================\n');
    assert(exist(mlPath,'file')==2, 'movieList not found: %s', mlPath);

    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);

    %% A) Segmentation (cell channel)
    fprintf('[%s] Step A: Segmentation | %d movies | chan=%d | mode=%s\n', ...
        datestr(now,'HH:MM:SS'), nMovies, tfmParams.segChan, tfmParams.segMode);
    drawnow;

    hSeg = waitbar(0, sprintf('Segmentation: 0/%d', nMovies), ...
        'Name', sprintf('ML %d/%d', ci, nML));

    for i = 1:nMovies
        MD = ML.getMovie(i);
        [~, mdName] = fileparts(MD.outputDirectory_);
        fprintf('[%s]   Seg %d/%d: %s\n', datestr(now,'HH:MM:SS'), i, nMovies, mdName);
        drawnow;

        tSeg = tic;
        runThresholdSegmentationForMD(MD, tfmParams);
        fprintf('[%s]   Seg %d/%d done (%.1f s)\n', ...
            datestr(now,'HH:MM:SS'), i, nMovies, toc(tSeg));

        if ishandle(hSeg)
            waitbar(i/nMovies, hSeg, ...
                sprintf('Segmentation: %d/%d  (%.0f s)', i, nMovies, toc(tSeg)));
        end
        drawnow;
    end
    if ishandle(hSeg), close(hSeg); end
    fprintf('[%s] Step A done.\n', datestr(now,'HH:MM:SS'));

    %% B) TFMPackage setup
    fprintf('[%s] Step B: TFMPackage setup (beadChan=%d)...\n', ...
        datestr(now,'HH:MM:SS'), tfmParams.iBeadChan);
    drawnow;
    setupTFMPackageForMovieList(mlPath, tfmParams);
    fprintf('[%s] Step B done.\n', datestr(now,'HH:MM:SS'));
    drawnow;

    %% C) Run TFM
    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);
    fprintf('[%s] Step C: tfmRun movie-by-movie | %d movies\n', ...
        datestr(now,'HH:MM:SS'), nMovies);
    fprintf('  This step runs displacement + force calculation.\n');
    fprintf('  Progress is printed per movie. Use Ctrl+C between movies if needed.\n');
    drawnow;

    hTFM = waitbar(0, 'TFM running... see Command Window for detail', ...
        'Name', sprintf('TFM ML %d/%d', ci, nML));
    drawnow;

    for i = 1:nMovies
        MD = ML.getMovie(i);
        [~, mdName] = fileparts(MD.outputDirectory_);
        fprintf('[%s] TFM %d/%d: %s\n', datestr(now,'HH:MM:SS'), i, nMovies, mdName);
        drawnow;
        tTFM = tic;
        try
            tfmRun(MD);
        catch ME
            warning(ME.identifier, 'tfmRun failed for %s: %s', mdName, ME.message);
        end
        elapsed = toc(tTFM);
        fprintf('[%s] TFM %d/%d done (%.1f s)\n', ...
            datestr(now,'HH:MM:SS'), i, nMovies, elapsed);
        if ishandle(hTFM)
            waitbar(i/nMovies, hTFM, ...
                sprintf('TFM: %d/%d  %s  (%.0f s)', i, nMovies, mdName, elapsed));
        end
        drawnow;
    end
    if ishandle(hTFM), close(hTFM); end

    fprintf('[%s] MovieList %d/%d DONE (total %.1f min)\n', ...
        datestr(now,'HH:MM:SS'), ci, nML, toc(tML)/60);
end

fprintf('\n[%s] All %d MovieLists finished.\n', datestr(now,'HH:MM:SS'), nML);

%% =====================================================================
%% LOCAL FUNCTIONS
%% =====================================================================

function MLpaths = pickMovieListsSafeUigetfile()
% Old-style file picker. This intentionally avoids uifigure/uilistbox.
MLpaths = {};
[fn, fp] = uigetfile( ...
    {'movieList.mat','movieList.mat'; '*.mat','MAT files (*.mat)'}, ...
    'Select movieList.mat file(s)', 'MultiSelect', 'on');
if isequal(fn,0)
    return;
end
if ischar(fn)
    fn = {fn};
end
MLpaths = cellfun(@(n) fullfile(fp,n), fn, 'UniformOutput', false);
MLpaths = unique(MLpaths(:)', 'stable');

bad = cellfun(@(p) exist(p,'file')~=2, MLpaths);
if any(bad)
    warning('run_TFM:missingML', 'Removing non-existent paths:\n%s', ...
        strjoin(MLpaths(bad), newline));
    MLpaths = MLpaths(~bad);
end
end

% ------------------------------------------------------------------
function ensureSegmentationPackageOnPath()
if exist('SegmentationPackage','class') == 8, return; end
warning('run_TFM:noSegPkg', 'SegmentationPackage not found on path. Please select its folder.');
root = uigetdir(pwd, 'Select folder containing SegmentationPackage');
assert(~isequal(root,0), 'User cancelled SegmentationPackage folder selection.');
addpath(genpath(root));
assert(exist('SegmentationPackage','class') == 8, ...
    'SegmentationPackage still not found after adding: %s', root);
fprintf('[Path] SegmentationPackage added: %s\n', root);
end

% ------------------------------------------------------------------
function runThresholdSegmentationForMD(MD, tfmParams)
% Configure and run SegmentationPackage on the CELL channel (segChan).
%   Process 1: ThresholdProcess  - automatic thresholding
%   Process 2: MaskRefinementProcess - morphological cleanup (default params)
%
% beadChan is NOT involved here.
%
% Supported ThresholdProcess modes (all automatic):
%   MinMax   : MethodIndx=1 - uses image min/max for threshold estimate
%   Otsu     : MethodIndx=2 - Otsu's method (default)
%   Rosin    : MethodIndx=3 - Rosin unimodal thresholding
%   Gradient : MethodIndx=4 - gradient-based thresholding

if ~tfmParams.forceRerunSeg
    anyMask = ~isempty(dir(fullfile(MD.outputDirectory_, 'masks', '**', '*.tif')));
    if anyMask
        fprintf('  [Seg] Skip (masks exist): %s\n', MD.outputDirectory_);
        return;
    end
end

iSeg = MD.getPackageIndex('SegmentationPackage');
if isempty(iSeg)
    MD.addPackage(SegmentationPackage(MD));
    iSeg = MD.getPackageIndex('SegmentationPackage');
end
segPack = MD.getPackage(iSeg);

% Process 1: ThresholdProcess
if isempty(segPack.processes_{1})
    segPack.createDefaultProcess(1);
end
proc1 = segPack.getProcess(1);
fp1   = proc1.funParams_;

fp1.ChannelIndex = tfmParams.segChan;
fp1.BatchMode    = true;

if isfield(fp1,'GaussFilterSigma')
    fp1.GaussFilterSigma = tfmParams.segGaussSigma;
end

segMode = lower(string(tfmParams.segMode));
switch segMode
    case "minmax"
        fp1.MethodIndx     = 1;
        fp1.PreThreshold   = false;
        fp1.IsPercentile   = false;
        fp1.ThresholdValue = [];
    case "otsu"
        fp1.MethodIndx     = 2;
        fp1.PreThreshold   = false;
        fp1.IsPercentile   = false;
        fp1.ThresholdValue = [];
    case "rosin"
        fp1.MethodIndx     = 3;
        fp1.PreThreshold   = false;
        fp1.IsPercentile   = false;
        fp1.ThresholdValue = [];
    case "gradient"
        fp1.MethodIndx     = 4;
        fp1.PreThreshold   = false;
        fp1.IsPercentile   = false;
        fp1.ThresholdValue = [];
    otherwise
        warning('run_TFM:unknownSegMode', ...
            'Unknown segMode "%s". Falling back to Otsu.', segMode);
        fp1.MethodIndx     = 2;
        fp1.PreThreshold   = false;
        fp1.IsPercentile   = false;
        fp1.ThresholdValue = [];
end

proc1.setPara(fp1);

fprintf('  [Seg-1/ThresholdProcess] chan=%d | mode=%s | sigma=%.1f | %s\n', ...
    tfmParams.segChan, tfmParams.segMode, tfmParams.segGaussSigma, MD.outputDirectory_);

if ismethod(proc1,'run')
    proc1.run();
else
    feval(proc1.funName_, MD, proc1.funParams_);
end

% Process 2: MaskRefinementProcess
try
    if numel(segPack.processes_) < 2 || isempty(segPack.processes_{2})
        segPack.createDefaultProcess(2);
    end
    proc2 = segPack.getProcess(2);
    fp2   = proc2.funParams_;

    if isfield(fp2,'ChannelIndex')
        fp2.ChannelIndex = tfmParams.segChan;
    end
    if isfield(fp2,'BatchMode')
        fp2.BatchMode = true;
    end

    proc2.setPara(fp2);

    fprintf('  [Seg-2/MaskRefinementProcess] chan=%d | default params | %s\n', ...
        tfmParams.segChan, MD.outputDirectory_);

    if ismethod(proc2,'run')
        proc2.run();
    else
        feval(proc2.funName_, MD, proc2.funParams_);
    end
catch ME
    warning(ME.identifier, ...
        'MaskRefinementProcess failed (non-fatal, ThresholdProcess mask still valid): %s', ME.message);
end

MD.save();
end

% ------------------------------------------------------------------
function tfmParams = promptTFMSegParams(defaults)
tfmParams = defaults;
fprintf('\n=== TFM Pipeline Settings ===\n');
fprintf('(Press Enter to keep default)\n\n');

fprintf('-- TFM (bead channel) --\n');
tmp = input(sprintf('Bead channel index [default %d]: ', defaults.iBeadChan), 's');
if ~isempty(tmp), tfmParams.iBeadChan = str2double(tmp); end

tmp = input(sprintf('Young''s modulus in Pa [default %.3g]: ', defaults.YoungModulus), 's');
if ~isempty(tmp), tfmParams.YoungModulus = str2double(tmp); end

tmp = input(sprintf('Poisson ratio [default %.2f]: ', defaults.PoissonRatio), 's');
if ~isempty(tmp), tfmParams.PoissonRatio = str2double(tmp); end

tmp = input(sprintf('Regularization param [default %.3g]: ', defaults.regParam), 's');
if ~isempty(tmp), tfmParams.regParam = str2double(tmp); end

fprintf('\n-- Displacement field tracking (proc2) --\n');
tmp = input(sprintf('Template window half-width in pixels [default %d]: ', defaults.minCorLength), 's');
if ~isempty(tmp), tfmParams.minCorLength = str2double(tmp); end

tmp = input(sprintf('Max displacement per frame in pixels [default %d]: ', defaults.maxFlowSpeed), 's');
if ~isempty(tmp), tfmParams.maxFlowSpeed = str2double(tmp); end

tmp = input(sprintf('Add non-local-max beads? 0/1 [default %d]: ', defaults.addNonLocMaxBeads), 's');
if ~isempty(tmp), tfmParams.addNonLocMaxBeads = logical(str2double(tmp)); end

fprintf('\n-- Segmentation (cell channel) --\n');
tmp = input(sprintf('Cell channel index for segmentation [default %d]: ', defaults.segChan), 's');
if ~isempty(tmp), tfmParams.segChan = str2double(tmp); end

modes = {'MinMax','Otsu','Rosin','Gradient'};
fprintf('Segmentation mode (all automatic):\n');
fprintf('  1) MinMax    - uses image min/max\n');
fprintf('  2) Otsu      - Otsu method (default)\n');
fprintf('  3) Rosin     - unimodal thresholding\n');
fprintf('  4) Gradient  - gradient-based\n');
tmp = input(sprintf('Choose [1-4] [default %d]: ', defaults.segModeIdx), 's');
if isempty(tmp)
    tfmParams.segMode = modes{defaults.segModeIdx};
else
    idx = max(1, min(4, round(str2double(tmp))));
    tfmParams.segMode = modes{idx};
end

tmp = input(sprintf('Gaussian blur sigma in pixels (0=none) [default %.1f]: ', defaults.segGaussSigma), 's');
if ~isempty(tmp), tfmParams.segGaussSigma = str2double(tmp); end

tmp = input(sprintf('Force rerun segmentation? 0/1 [default %d]: ', defaults.forceRerunSeg), 's');
if ~isempty(tmp), tfmParams.forceRerunSeg = logical(str2double(tmp)); end

fprintf('\n=== Summary ===\n');
fprintf('  beadChan       = %d  (drift correction + displacement tracking)\n', tfmParams.iBeadChan);
fprintf('  segChan        = %d  (ThresholdProcess + MaskRefinementProcess)\n', tfmParams.segChan);
fprintf('  E              = %.3g Pa | nu=%.2f | regParam=%.3g\n', ...
    tfmParams.YoungModulus, tfmParams.PoissonRatio, tfmParams.regParam);
fprintf('  minCorLength   = %d px (template window)\n', tfmParams.minCorLength);
fprintf('  maxFlowSpeed   = %d px/frame\n',              tfmParams.maxFlowSpeed);
fprintf('  addNonLocMax   = %d\n',                       tfmParams.addNonLocMaxBeads);
fprintf('  segMode        = %s | sigma=%.1f px | forceRerun=%d\n\n', ...
    tfmParams.segMode, tfmParams.segGaussSigma, tfmParams.forceRerunSeg);
end
