%% run_TFMpipeline_interactive.m
% Segmentation (Threshold) -> TFM pipeline runner (interactive, reproducible)
%
% What it does
%   1) Lets you select one or more MovieList files (movieList.mat)
%   2) Runs cell segmentation FIRST using SegmentationPackage / ThresholdProcess
%      - Choose segmentation mode (Manual / Otsu / Adaptive-like)
%      - Choose segmentation channel
%      - Optional Gaussian blur (sigma in pixels) prior to thresholding
%   3) Sets up TFMPackage using your existing setupTFMPackageForMovieList.m
%      - Uses interactive bead channel, Young's modulus, Poisson ratio, regParam
%   4) Ensures Step-4 (ForceFieldCalculationProcess) regParam is set to the chosen value
%   5) Runs tfmRunML (serial by default)
%
% Notes
%   - This script tries to use real ThresholdProcess funParams_ fields safely:
%       ChannelIndex, MethodIndx, ThresholdValue, IsPercentile, PreThreshold,
%       GaussFilterSigma, BatchMode.
%   - If some fields differ in your local ThresholdProcess implementation,
%     the script will fall back gracefully where possible.
%
% Sangyoon Han (cleaned/standardized version)

clear; clc;

%% ===== USER SETTINGS =====
useParallel = false;     % default: serial
maxWorkers  = 12;        % only used if useParallel=true

%% ===== DEFAULTS (interactive prompts) =====
defaults = struct();

% ---- TFM defaults ----
defaults.iBeadChan    = 3;
defaults.YoungModulus = 2700;   % Pa
defaults.PoissonRatio = 0.49;
defaults.regParam     = 0.005;

% ---- Segmentation defaults ----
defaults.segModeIdx     = 2;    % 1 ManualThreshold, 2 Otsu, 3 Adaptive
defaults.manualThresh   = 0.0;  % for "whole FOV is cell" case
defaults.segChan        = 1;
defaults.segGaussSigma  = 2;    % pixels (0 = none). You typically use 2~3.

% Prompt user for parameters
tfmParams = promptTFMParams(defaults);

% Derive segmentation inputs
segChan       = tfmParams.segChan;
segInvert     = false;  % keep for compatibility; most ThresholdProcess versions don't support it directly
segForceRerun = false;  % if masks already exist, skip unless true

%% ===== PICK MOVIELISTS =====
MLpaths = pickMovieListsMultiFolderGUI();
assert(~isempty(MLpaths), 'No movieList.mat selected.');

%% ===== Ensure SegmentationPackage on path (ask once if missing) =====
ensureSegmentationPackageOnPath();

%% ===== Optional: parallel pool =====
useParallel = ensureParallelPool(useParallel, maxWorkers);

%% ===== MAIN LOOP =====
for ci = 1:numel(MLpaths)

    mlPath = char(MLpaths{ci});

    fprintf('\n====================================================\n');
    fprintf('MovieList %d/%d\nML: %s\n', ci, numel(MLpaths), mlPath);
    fprintf('====================================================\n');
    assert(exist(mlPath,'file')==2, 'movieList not found: %s', mlPath);

    %% A) Load ML and run segmentation first
    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);
    fprintf('[Info] nMovies = %d\n', nMovies);

    fprintf('[Seg] Running Threshold segmentation (chan=%d, mode=%s, gauss=%.2f px)...\n', ...
        segChan, tfmParams.segMode, tfmParams.segGaussSigma);

    for i = 1:nMovies
        MD = ML.getMovie(i);
        runThresholdSegmentationForMD(MD, segChan, segInvert, segForceRerun, tfmParams);
    end
    fprintf('[Seg] Done.\n');

    %% B) Setup TFMPackage (uses your existing setup code)
    fprintf('[TFM-Setup] Setting up TFMPackage...\n');
    setupTFMPackageForMovieList(mlPath, tfmParams);

    %% C) Reload ML after setup
    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);

    %% D) Force Step-4 regParam to chosen value
    fprintf('[TFM] Setting Step-4 regParam = %.6g\n', tfmParams.regParam);

    for i = 1:nMovies
        MD = ML.getMovie(i);

        iPack = MD.getPackageIndex('TFMPackage');
        assert(~isempty(iPack), 'TFMPackage not found in MD: %s', MD.outputDirectory_);

        pack = MD.getPackage(iPack);
        assert(numel(pack.processes_) >= 4 && ~isempty(pack.processes_{4}), ...
            'TFM Step-4 process missing in MD: %s', MD.outputDirectory_);

        forceProc = pack.processes_{4};
        p = forceProc.funParams_;

        if isfield(p,'regParam')
            p.regParam = tfmParams.regParam;
        elseif isfield(p,'regularizationParameter')
            p.regularizationParameter = tfmParams.regParam;
        else
            warning('No regParam field found in Step-4 funParams_ for %s', MD.outputDirectory_);
        end

        forceProc.setPara(p);
        MD.save();
    end

    %% E) Run TFM
    fprintf('[TFM] Running tfmRunML (parallel=%d)...\n', useParallel);
    tfmRunML(ML, useParallel);

    fprintf('[Done] %s\n', mlPath);
end

disp('All MovieLists finished.');

%% ===================== Segmentation helpers =====================

function ensureSegmentationPackageOnPath()
% Checks if SegmentationPackage is available; if not, prompts for folder and adds it.

if exist('SegmentationPackage','class') == 8
    return;
end

warning(['SegmentationPackage not found on path.' newline ...
         'Please select the folder that contains the SegmentationPackage code.']);
root = uigetdir(pwd, 'Select folder containing SegmentationPackage (extracted zip folder)');
assert(~isequal(root,0), 'User canceled selecting segmentation package folder.');

addpath(genpath(root));

assert(exist('SegmentationPackage','class') == 8, ...
    'Still cannot find SegmentationPackage after addpath(genpath(%s)). Check folder selection.', root);

fprintf('[Path] Added segmentation package path: %s\n', root);
end

function runThresholdSegmentationForMD(MD, segChan, invert, forceRerun, tfmParams)
% Adds SegmentationPackage + ThresholdProcess if needed, configures it using
% ThresholdProcess funParams_ fields (robustly), and runs it.

% ---- skip if masks exist ----
if ~forceRerun
    anyMask = ~isempty(dir(fullfile(MD.outputDirectory_, 'masks', '**', '*.tif')));
    if anyMask
        fprintf('  [Seg] Skip (masks already exist): %s\n', MD.outputDirectory_);
        return;
    end
end

% ---- Ensure SegmentationPackage ----
iSeg = MD.getPackageIndex('SegmentationPackage');
if isempty(iSeg)
    MD.addPackage(SegmentationPackage(MD));
    iSeg = MD.getPackageIndex('SegmentationPackage');
end
segPack = MD.getPackage(iSeg);

% ---- Ensure ThresholdProcess (process #1 in SegmentationPackage) ----
if isempty(segPack.processes_{1})
    segPack.createDefaultProcess(1);
end
proc = segPack.getProcess(1);

% ---- Get funParams ----
fp = proc.funParams_;

% Segment only segChan (scalar)
if isfield(fp,'ChannelIndex')
    fp.ChannelIndex = segChan;
end

% Make it non-interactive (batch)
if isfield(fp,'BatchMode')
    fp.BatchMode = true;
end

% Gaussian blur before thresholding (sigma in pixels)
if isfield(fp,'GaussFilterSigma')
    fp.GaussFilterSigma = tfmParams.segGaussSigma;
end

% Invert: many ThresholdProcess versions do not support it; warn once per MD
if invert
    warning('[Seg] invert=true requested, but no standard ThresholdProcess parameter is enforced here. Ignoring invert.');
end

% ---- Configure method based on tfmParams.segMode ----
segMode = lower(string(tfmParams.segMode));

% Many HanLab ThresholdProcess implementations use:
%   MethodIndx: 1=MinMax, 2=Otsu, 3=Rosin, 4=Gradient-based
% We set those when available; otherwise, we try "Method" string if present.
if isfield(fp,'MethodIndx')
    switch segMode
        case "manualthreshold"
            fp.MethodIndx = 1;  % MinMax (but using ThresholdValue you supply)
        case "otsu"
            fp.MethodIndx = 2;
        case "adaptive"
            fp.MethodIndx = 4;  % Gradient-based as the closest "adaptive-like" option
        otherwise
            fp.MethodIndx = 2;
    end
elseif isfield(fp,'Method')
    switch segMode
        case "manualthreshold"
            fp.Method = 'Manual';
        case "otsu"
            fp.Method = 'Otsu';
        case "adaptive"
            fp.Method = 'Gradient';
        otherwise
            fp.Method = 'Otsu';
    end
end

% Manual threshold fields (only meaningful in manual mode)
if segMode == "manualthreshold"
    if isfield(fp,'ThresholdValue')
        fp.ThresholdValue = tfmParams.manualThresh;
    end
    if isfield(fp,'IsPercentile')
        fp.IsPercentile = false;
    end
    if isfield(fp,'PreThreshold')
        fp.PreThreshold = false;
    end
else
    % Let the method decide threshold unless your implementation expects empty vs NaN.
    if isfield(fp,'ThresholdValue')
        fp.ThresholdValue = [];  %#ok<NASGU>
    end
end

% ---- Apply updated params ----
proc.setPara(fp);

% ---- Run ----
fprintf('  [Seg] Run ThresholdProcess: %s (chan=%d, mode=%s, gauss=%.2f)\n', ...
    MD.outputDirectory_, segChan, tfmParams.segMode, tfmParams.segGaussSigma);

if ismethod(proc,'run')
    proc.run();
else
    feval(proc.funName_, MD, proc.funParams_);
end

MD.save();
end

%% ===================== Helper: MovieList picker =====================

function MLpaths = pickMovieListsMultiFolderGUI()
% Multi-folder MovieList picker (stable shape / duplicates handled)

if ~usejava('desktop')
    error('GUI not available (no desktop). Use CLI entry instead.');
end

MLpaths = {};

f = uifigure('Name','Select MovieLists (multi-folder)','Position',[200 200 860 420]);

lst = uilistbox(f, ...
    'Position',[20 70 820 330], ...
    'Items',{}, ...
    'Multiselect','on');

uibutton(f,'Text','Add MovieList(s)...','Position',[20 20 140 30], ...
    'ButtonPushedFcn', @(~,~) addFiles());

uibutton(f,'Text','Add from Folder (recursive)...','Position',[170 20 190 30], ...
    'ButtonPushedFcn', @(~,~) addFromFolderRecursive());

uibutton(f,'Text','Remove selected','Position',[370 20 130 30], ...
    'ButtonPushedFcn', @(~,~) removeSelected());

uibutton(f,'Text','Clear','Position',[510 20 80 30], ...
    'ButtonPushedFcn', @(~,~) clearAll());

uibutton(f,'Text','OK','Position',[720 20 50 30], ...
    'ButtonPushedFcn', @(~,~) uiresume(f));

uibutton(f,'Text','Cancel','Position',[780 20 60 30], ...
    'ButtonPushedFcn', @(~,~) cancel());

uiwait(f);

if ~isvalid(f)
    MLpaths = {};
    return;
end

items = string(lst.Items(:));
delete(f);

MLpaths = uniqueStableCell(cellstr(items(:)));

% validate paths
bad = cellfun(@(p) exist(p,'file')~=2, MLpaths);
if any(bad)
    warning('Some selected paths do not exist and will be removed:\n%s', strjoin(string(MLpaths(bad)), newline));
    MLpaths = MLpaths(~bad);
end

    function addFiles()
        [fn, fp] = uigetfile({'movieList.mat','movieList.mat'; '*.mat','MAT-files (*.mat)'}, ...
            'Select movieList.mat file(s)', 'MultiSelect','on');
        if isequal(fn,0), return; end
        if ischar(fn) || isstring(fn), fn = {fn}; end

        newPaths = strings(numel(fn),1);
        for ii = 1:numel(fn)
            newPaths(ii) = fullfile(fp, fn{ii});
        end

        curr = string(lst.Items(:));
        merged = uniqueStableStr([curr(:); newPaths(:)]);
        lst.Items = cellstr(merged(:));
    end

    function addFromFolderRecursive()
        root = uigetdir(pwd, 'Pick a folder (search recursively for movieList.mat)');
        if isequal(root,0), return; end

        hits = recursiveFindMovieLists(root);
        if isempty(hits)
            uialert(f, sprintf('No movieList.mat found under:\n%s', root), 'No matches');
            return;
        end

        curr = string(lst.Items(:));
        merged = uniqueStableStr([curr(:); string(hits(:))]);
        lst.Items = cellstr(merged(:));
    end

    function removeSelected()
        sel = string(lst.Value(:));
        if isempty(sel), return; end
        curr = string(lst.Items(:));
        keep = ~ismember(curr, sel);
        lst.Items = cellstr(curr(keep));
    end

    function clearAll()
        lst.Items = {};
    end

    function cancel()
        delete(f);
        MLpaths = {};
    end
end

function hits = recursiveFindMovieLists(rootFolder)
d = dir(fullfile(rootFolder, '**', 'movieList.mat'));
d = d(~[d.isdir]);
hits = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
end

function out = uniqueStableStr(s)
s = string(s(:));
[~, ia] = unique(s, 'stable');
out = s(ia);
end

function out = uniqueStableCell(c)
s = string(c(:));
[~, ia] = unique(s, 'stable');
out = cellstr(s(ia));
end

%% ===================== Interactive prompts =====================

function tfmParams = promptTFMParams(defaults)
% Interactive prompts (works in terminal MATLAB too)

tfmParams = defaults;

fprintf('\n=== Interactive settings ===\n');

% ---- bead channel ----
tmp = input(sprintf('Bead channel index (default %d): ', defaults.iBeadChan), 's');
if ~isempty(tmp), tfmParams.iBeadChan = str2double(tmp); end

% ---- Young modulus ----
tmp = input(sprintf('Young''s modulus in Pa (default %.6g): ', defaults.YoungModulus), 's');
if ~isempty(tmp), tfmParams.YoungModulus = str2double(tmp); end

% ---- Poisson ratio ----
tmp = input(sprintf('Poisson ratio (default %.2f): ', defaults.PoissonRatio), 's');
if ~isempty(tmp), tfmParams.PoissonRatio = str2double(tmp); end

% ---- regularization parameter ----
tmp = input(sprintf('TFM regParam (default %.6g): ', defaults.regParam), 's');
if ~isempty(tmp), tfmParams.regParam = str2double(tmp); end

% ---- segmentation mode ----
modes = {'ManualThreshold','Otsu','Adaptive'};
fprintf('\nSegmentation mode options:\n');
for i=1:numel(modes)
    fprintf('  %d) %s\n', i, modes{i});
end
tmp = input(sprintf('Choose segmentation mode [1-%d] (default %d): ', numel(modes), defaults.segModeIdx), 's');
if isempty(tmp)
    tfmParams.segMode = modes{defaults.segModeIdx};
else
    idx = max(1, min(numel(modes), round(str2double(tmp))));
    tfmParams.segMode = modes{idx};
end

% ---- manual threshold if chosen ----
if strcmpi(tfmParams.segMode,'ManualThreshold')
    tmp = input(sprintf('Manual mask threshold (0~1) (default %.3f): ', defaults.manualThresh), 's');
    if ~isempty(tmp), tfmParams.manualThresh = str2double(tmp); end
end

% ---- segmentation channel ----
tmp = input(sprintf('Segmentation channel index (default %d): ', defaults.segChan), 's');
if ~isempty(tmp), tfmParams.segChan = str2double(tmp); end

% ---- Gaussian blur sigma ----
tmp = input(sprintf('Segmentation Gaussian sigma in pixels (0 = none, default %.2f): ', defaults.segGaussSigma), 's');
if ~isempty(tmp)
    tfmParams.segGaussSigma = str2double(tmp);
else
    tfmParams.segGaussSigma = defaults.segGaussSigma;
end

fprintf('\n[TFM] beadChan=%d | E=%.6g Pa | nu=%.2f | reg=%.6g\n', ...
    tfmParams.iBeadChan, tfmParams.YoungModulus, tfmParams.PoissonRatio, tfmParams.regParam);
fprintf('[Seg] mode=%s | chan=%d | gauss=%.2f px\n\n', ...
    tfmParams.segMode, tfmParams.segChan, tfmParams.segGaussSigma);
end

%% ===================== Parallel helper =====================

function useParallel = ensureParallelPool(useParallel, maxWorkers)
if ~useParallel
    return;
end

try
    if license('test','Distrib_Computing_Toolbox')
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local', maxWorkers);
        end
    else
        warning('Parallel Computing Toolbox not available. Using serial.');
        useParallel = false;
    end
catch ME
    warning(ME.identifier,'Parallel pool failed (%s). Using serial.', ME.message);
    useParallel = false;
end
end