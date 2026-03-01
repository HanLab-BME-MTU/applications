%% run_TFMpipeline_interactive.m
% Segmentation (ThresholdProcess) -> TFMPackage pipeline runner
%   - Select one or more movieList.mat files (multi-folder GUI)
%   - Run segmentation first (SegmentationPackage: ThresholdProcess)
%       * Supports Manual / Otsu / Gradient-based
%       * Optional Gaussian blur (sigma in pixels)
%   - Setup TFMPackage using user-provided parameters (no hardcoding)
%   - Optionally set Step-4 regularization parameter
%   - Run tfmRunML (serial by default)
%
% Sangyoon Han (cleaned/updated)

clear; clc;

%% ===== USER SETTINGS =====
useParallel = false;   % recommended: serial for reproducibility

% ---- Defaults (user can override interactively) ----
defaults = struct();
defaults.iBeadChan       = 3;
defaults.YoungModulus    = 2700;   % Pa
defaults.PoissonRatio    = 0.49;
defaults.regParam        = 0.005;

defaults.segChan         = 1;      % channel used for segmentation
defaults.segModeIdx      = 2;      % 1 Manual, 2 Otsu, 3 Gradient-based
defaults.manualThresh    = 0.0;    % used only if Manual
defaults.segGaussSigma   = 2;      % pixels (0 = none)

defaults.forceRerunSeg   = false;  % if masks exist, skip unless true

tfmParams = promptTFMSegParams(defaults);

%% ===== PICK MOVIELISTS =====
MLpaths = pickMovieListsMultiFolderGUI();
assert(~isempty(MLpaths), 'No movieList.mat selected.');

%% ===== Ensure SegmentationPackage on path (ask once if missing) =====
ensureSegmentationPackageOnPath();

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

    fprintf('[Seg] Running ThresholdProcess (chan=%d, mode=%s, GaussSigma=%.2f px)\n', ...
        tfmParams.segChan, tfmParams.segMode, tfmParams.segGaussSigma);

    for i = 1:nMovies
        MD = ML.getMovie(i);
        runThresholdSegmentationForMD(MD, tfmParams);
        progressText(i/nMovies, sprintf('Segmentation: %d/%d', i, nMovies));
    end
    fprintf('\n[Seg] Done.\n');

    %% B) Setup TFMPackage
    fprintf('[TFM-Setup] Setting up TFMPackage...\n');
    setupTFMPackageForMovieList(mlPath, tfmParams);

    %% C) Reload ML after setup
    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);

    %% D) Set Step-4 regParam (optional, but usually desired)
    fprintf('[TFM] Setting Step-4 regParam = %.6g\n', tfmParams.regParam);
    for i = 1:nMovies
        MD = ML.getMovie(i);
        iPack = MD.getPackageIndex('TFMPackage');
        pack = MD.getPackage(iPack);

        if numel(pack.processes_) >= 4 && ~isempty(pack.processes_{4})
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
        else
            warning('TFMPackage Step-4 not found for %s', MD.outputDirectory_);
        end
    end

    %% E) Run TFM
    fprintf('[TFM] Running tfmRunML (parallel=%d)...\n', useParallel);
    tfmRunML(ML, useParallel);

    fprintf('[Done] %s\n', mlPath);
end

disp('All MovieLists finished.');

%% ===================== Segmentation helpers =====================

function ensureSegmentationPackageOnPath()
% Ensure SegmentationPackage is available on MATLAB path.

if exist('SegmentationPackage','class') == 8
    return;
end

warning('SegmentationPackage not found on path. Please select the folder that contains SegmentationPackage code.');
root = uigetdir(pwd, 'Select folder containing SegmentationPackage');
assert(~isequal(root,0), 'User canceled selecting segmentation package folder.');

addpath(genpath(root));

assert(exist('SegmentationPackage','class') == 8, ...
    'Still cannot find SegmentationPackage after addpath(genpath(%s)). Check folder selection.', root);

fprintf('[Path] Added SegmentationPackage path: %s\n', root);
end

function runThresholdSegmentationForMD(MD, tfmParams)
% Configure and run ThresholdProcess using REAL ThresholdProcess fields.

% ---- Skip if masks exist ----
if ~tfmParams.forceRerunSeg
    anyMask = ~isempty(dir(fullfile(MD.outputDirectory_, 'masks', '**', '*.tif')));
    if anyMask
        fprintf('  [Seg] Skip (masks exist): %s\n', MD.outputDirectory_);
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

% ---- Ensure ThresholdProcess (process #1) ----
if isempty(segPack.processes_{1})
    segPack.createDefaultProcess(1);
end
proc = segPack.getProcess(1);

% ---- Get & set funParams ----
fp = proc.funParams_;

% Always segment only one channel here
fp.ChannelIndex = tfmParams.segChan;

% No GUI / batch mode
fp.BatchMode = true;

% Optional Gaussian blur sigma in pixels (0 = none)
if isfield(fp,'GaussFilterSigma')
    fp.GaussFilterSigma = tfmParams.segGaussSigma;
end

% Select method
% ThresholdProcess supports:
%   1 = MinMax (manual threshold)
%   2 = Otsu
%   3 = Rosin
%   4 = Gradient (best ?adaptive-ish? option in this package)
segMode = lower(string(tfmParams.segMode));

switch segMode
    case "manual"
        fp.MethodIndx      = 1;
        fp.PreThreshold    = false;
        fp.IsPercentile    = false;
        fp.ThresholdValue  = tfmParams.manualThresh;  % e.g., 0 for "whole FOV is cell"

    case "otsu"
        fp.MethodIndx      = 2;
        fp.PreThreshold    = false;
        fp.IsPercentile    = false;
        fp.ThresholdValue  = [];   % let method decide

    case "gradient"
        fp.MethodIndx      = 4;
        fp.PreThreshold    = false;
        fp.IsPercentile    = false;
        fp.ThresholdValue  = [];   % let method decide

    otherwise
        warning('[Seg] Unknown segMode=%s. Falling back to Otsu.', segMode);
        fp.MethodIndx      = 2;
        fp.PreThreshold    = false;
        fp.IsPercentile    = false;
        fp.ThresholdValue  = [];
end

proc.setPara(fp);

% ---- Run ----
fprintf('  [Seg] ThresholdProcess: %s (chan=%d, mode=%s, sigma=%.2f)\n', ...
    MD.outputDirectory_, tfmParams.segChan, tfmParams.segMode, tfmParams.segGaussSigma);

if ismethod(proc,'run')
    proc.run();
else
    feval(proc.funName_, MD, proc.funParams_);
end

MD.save();
end

%% ===================== MovieList picker (multi-folder GUI) =====================

function MLpaths = pickMovieListsMultiFolderGUI()
% Multi-folder MovieList picker.

if ~usejava('desktop')
    error('GUI not available (no desktop). Please select movieList.mat paths in CLI mode.');
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
    'ButtonPushedFcn', @(~,~) set(lst,'Items',{}));

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

MLpaths = cellstr(unique(items,'stable'));

% validate existence
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
        merged = unique([curr; newPaths(:)], 'stable');
        lst.Items = cellstr(merged(:));
    end

    function addFromFolderRecursive()
        root = uigetdir(pwd, 'Pick a folder (search recursively for movieList.mat)');
        if isequal(root,0), return; end

        d = dir(fullfile(root, '**', 'movieList.mat'));
        hits = arrayfun(@(x) fullfile(x.folder, x.name), d(~[d.isdir]), 'UniformOutput', false);
        if isempty(hits)
            uialert(f, sprintf('No movieList.mat found under:\n%s', root), 'No matches');
            return;
        end

        curr = string(lst.Items(:));
        merged = unique([curr; string(hits(:))], 'stable');
        lst.Items = cellstr(merged(:));
    end

    function removeSelected()
        sel = string(lst.Value(:));
        if isempty(sel), return; end
        curr = string(lst.Items(:));
        lst.Items = cellstr(curr(~ismember(curr, sel)));
    end

    function cancel()
        delete(f);
        MLpaths = {};
    end
end

%% ===================== Interactive prompt for params =====================

function tfmParams = promptTFMSegParams(defaults)

tfmParams = defaults;

fprintf('\n=== Interactive settings ===\n');

% Bead channel
tmp = input(sprintf('Bead channel index (default %d): ', defaults.iBeadChan), 's');
if ~isempty(tmp), tfmParams.iBeadChan = str2double(tmp); end

% Young modulus
tmp = input(sprintf('Young''s modulus in Pa (default %.3g): ', defaults.YoungModulus), 's');
if ~isempty(tmp), tfmParams.YoungModulus = str2double(tmp); end

% Poisson ratio
tmp = input(sprintf('Poisson ratio (default %.2f): ', defaults.PoissonRatio), 's');
if ~isempty(tmp), tfmParams.PoissonRatio = str2double(tmp); end

% Regularization parameter
tmp = input(sprintf('TFM regParam (default %.3g): ', defaults.regParam), 's');
if ~isempty(tmp), tfmParams.regParam = str2double(tmp); end

% Segmentation channel
tmp = input(sprintf('Segmentation channel index (default %d): ', defaults.segChan), 's');
if ~isempty(tmp), tfmParams.segChan = str2double(tmp); end

% Segmentation mode (ThresholdProcess MethodIndx)
modes = {'Manual','Otsu','Gradient'};
fprintf('\nSegmentation mode options (ThresholdProcess):\n');
fprintf('  1) Manual (MinMax)  -> uses ThresholdValue\n');
fprintf('  2) Otsu             -> automatic\n');
fprintf('  3) Gradient-based   -> more adaptive-like\n');
tmp = input(sprintf('Choose segmentation mode [1-%d] (default %d): ', numel(modes), defaults.segModeIdx), 's');

if isempty(tmp)
    tfmParams.segMode = modes{defaults.segModeIdx};
else
    idx = max(1, min(numel(modes), round(str2double(tmp))));
    tfmParams.segMode = modes{idx};
end

% Manual threshold value if Manual selected
if strcmpi(tfmParams.segMode,'Manual')
    tmp = input(sprintf('Manual threshold value (default %.3f): ', defaults.manualThresh), 's');
    if ~isempty(tmp), tfmParams.manualThresh = str2double(tmp); end
else
    tfmParams.manualThresh = defaults.manualThresh;
end

% Gaussian blur sigma
tmp = input(sprintf('Gaussian blur sigma (pixels, 0 = none) (default %.2f): ', defaults.segGaussSigma), 's');
if ~isempty(tmp), tfmParams.segGaussSigma = str2double(tmp); end

% Force rerun segmentation
tmp = input(sprintf('Force rerun segmentation if masks exist? 0/1 (default %d): ', defaults.forceRerunSeg), 's');
if ~isempty(tmp), tfmParams.forceRerunSeg = logical(str2double(tmp)); end

fprintf('\n[TFM] beadChan=%d | E=%.3g Pa | nu=%.2f | reg=%.3g\n', ...
    tfmParams.iBeadChan, tfmParams.YoungModulus, tfmParams.PoissonRatio, tfmParams.regParam);
fprintf('[Seg] chan=%d | mode=%s | manualThr=%.3f | GaussSigma=%.2f px | forceRerun=%d\n\n', ...
    tfmParams.segChan, tfmParams.segMode, tfmParams.manualThresh, tfmParams.segGaussSigma, tfmParams.forceRerunSeg);

end