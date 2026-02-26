%% run_TFMpipeline_interactive.m
% Segmentation (Threshold) -> TFM pipeline runner
%   - pick one or more movieList.mat files (multi-folder GUI)
%   - run cell segmentation first (SegmentationPackage: ThresholdProcess)
%       * manual threshold = 0 (whole FOV cell)
%   - setup TFMPackage
%   - FORCE Step-4 regularization parameter = 0.005
%   - run tfmRunML (serial; until step 5)
%
% Sangyoon Han

clear; clc;

%% ===== USER SETTINGS =====
useParallel = false;        % serial (as requested)

% ---- TFM ----
forcedRegParam = 0.005;

tfmParams = struct();
tfmParams.beadChan     = 3;
tfmParams.YoungModulus = 2700;     % Pa
tfmParams.regParam     = forcedRegParam;

% ---- Segmentation (Threshold) ----
segChan = 1;                % <<< change if your "cell" channel is not 1
segThreshold = 0;           % <<< manual threshold
segInvert = false;          % usually false
segForceRerun = false;      % if masks already exist, skip unless true

%% ===== PICK MOVIELISTS =====
MLpaths = pickMovieListsMultiFolderGUI();
assert(~isempty(MLpaths), 'No movieList.mat selected.');

%% ===== Ensure SegmentationPackage on path (ask once if missing) =====
ensureSegmentationPackageOnPath();

%% ===== MAIN =====
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

    fprintf('[Seg] Running Threshold segmentation (chan=%d, thr=%g)...\n', segChan, segThreshold);
    for i = 1:nMovies
        MD = ML.getMovie(i);
        runThresholdSegmentationForMD(MD, segChan, segThreshold, segInvert, segForceRerun);
    end
    fprintf('[Seg] Done.\n');

    %% B) Setup TFMPackage
    fprintf('[TFM-Setup] Setting up TFMPackage...\n');
    setupTFMPackageForMovieList(mlPath, tfmParams);

    %% C) Reload ML after setup
    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);

    %% D) Force Step-4 regParam
    fprintf('[TFM] Forcing Step-4 regParam = %.4f\n', forcedRegParam);
    for i = 1:nMovies
        MD = ML.getMovie(i);
        iPack = MD.getPackageIndex('TFMPackage');
        pack = MD.getPackage(iPack);
        forceProc = pack.processes_{4};
        p = forceProc.funParams_;

        if isfield(p,'regParam')
            p.regParam = forcedRegParam;
        elseif isfield(p,'regularizationParameter')
            p.regularizationParameter = forcedRegParam;
        else
            warning('No regParam field found in Step-4 funParams_ for %s', MD.outputDirectory_);
        end

        forceProc.setPara(p);
        MD.save();
    end

    %% E) Run TFM (serial)
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

warning('SegmentationPackage not found on path. Please select the folder that contains "segmentationPackage" code.');
root = uigetdir(pwd, 'Select folder containing SegmentationPackage (the extracted zip folder)');
assert(~isequal(root,0), 'User canceled selecting segmentation package folder.');

addpath(genpath(root));

assert(exist('SegmentationPackage','class') == 8, ...
    'Still cannot find SegmentationPackage after addpath(genpath(%s)). Check folder selection.', root);

fprintf('[Path] Added segmentation package path: %s\n', root);
end

function runThresholdSegmentationForMD(MD, segChan, thr, invert, forceRerun)
% Adds SegmentationPackage + ThresholdProcess if needed, sets manual threshold, runs it.
% Writes masks under <MD.outputDirectory_>/masks (default behavior).

% If masks exist and not forcing rerun, skip
maskDir = fullfile(MD.outputDirectory_, 'masks', sprintf('channel_%d', segChan));
% NOTE: ThresholdProcess actually creates masks/channel_* or channel-specific dirs depending on implementation.
% We'll use a broader check for any tif in masks.
if ~forceRerun
    anyMask = ~isempty(dir(fullfile(MD.outputDirectory_, 'masks', '**', '*.tif')));
    if anyMask
        fprintf('  [Seg] Skip (masks already exist): %s\n', MD.outputDirectory_);
        return;
    end
end

% Ensure SegmentationPackage
iSeg = MD.getPackageIndex('SegmentationPackage');
if isempty(iSeg)
    MD.addPackage(SegmentationPackage(MD));
    iSeg = MD.getPackageIndex('SegmentationPackage');
end
segPack = MD.getPackage(iSeg);

% Ensure ThresholdProcess (process #1 in SegmentationPackage)
if isempty(segPack.processes_{1})
    segPack.createDefaultProcess(1);
end
proc = segPack.getProcess(1);

% Set parameters
fp = proc.funParams_;

% Segment only segChan
if isfield(fp,'ChannelIndex')
    fp.ChannelIndex = segChan;
end

% Manual thresholding
if isfield(fp,'ThresholdValue')
    fp.ThresholdValue = thr;
end
if isfield(fp,'IsPercentile')
    fp.IsPercentile = false;
end

% Make it non-interactive
if isfield(fp,'BatchMode')
    fp.BatchMode = true;
end

% Invert (rarely needed)
if isfield(fp,'Invert')
    fp.Invert = logical(invert);
end

% Ensure output directory exists
if isfield(fp,'OutputDirectory') && ~isempty(fp.OutputDirectory)
    if ~exist(fp.OutputDirectory,'dir'); mkdir(fp.OutputDirectory); end
end

proc.setPara(fp);

% Run
fprintf('  [Seg] Run ThresholdProcess: %s\n', MD.outputDirectory_);
if ismethod(proc,'run')
    proc.run();
else
    % Fallback: call process function directly
    feval(proc.funName_, MD, proc.funParams_);
end

MD.save();
end


%% ===================== Helper: pick movie lists =====================
function MLpaths = pickMovieListsInteractive()
% GUI multi-select for movieList.mat files.
% CLI fallback: prompt for N and paths.

MLpaths = {};

if usejava('desktop')
    [fn, fp] = uigetfile({'*.mat','MAT-files (*.mat)'}, ...
                         'Select one or more movieList.mat files', ...
                         'MultiSelect','on');
    if isequal(fn,0)
        error('Selection canceled.');
    end
    if ischar(fn) || isstring(fn)
        fn = {fn};
    end
    MLpaths = cell(numel(fn),1);
    for i=1:numel(fn)
        MLpaths{i} = fullfile(fp, fn{i});
    end

    % optional: filter only movieList.mat
    isML = endsWith(string(MLpaths), "movieList.mat", 'IgnoreCase', true);
    if any(~isML)
        warning('Some selected files are not named "movieList.mat". They will still be attempted.');
    end
    return;
end

% CLI fallback
fprintf('\n[Interactive] CLI mode: enter movieList.mat paths\n');
n = input('How many MovieLists? ');
if isempty(n) || ~isscalar(n) || n<1
    error('Invalid number.');
end
MLpaths = cell(n,1);
for i=1:n
    p = input(sprintf('Path %d/%d: ', i, n), 's');
    if isempty(p), error('Empty path.'); end
    MLpaths{i} = string(p);
end
end

function MLpaths = pickMovieListsMultiFolderGUI()
% Multi-folder MovieList picker (fixed shape handling)

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

items = string(lst.Items(:));   % column
delete(f);

MLpaths = cellstr(items);

% remove duplicates, keep order
MLpaths = uniqueStableCell(MLpaths);

% validate
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

        curr = string(lst.Items(:));      % force column
        newp = string(newPaths(:));       % force column
        merged = uniqueStableStr([curr; newp]);
        lst.Items = cellstr(merged(:));   % assign as cell vector
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
        newp = string(hits(:));
        merged = uniqueStableStr([curr; newp]);
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

%% ===== helper: recursive finder =====
function hits = recursiveFindMovieLists(rootFolder)
d = dir(fullfile(rootFolder, '**', 'movieList.mat'));
d = d(~[d.isdir]);
hits = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
end

%% ===== helper: stable unique (string array) =====
function out = uniqueStableStr(s)
s = string(s(:));
[~, ia] = unique(s, 'stable');
out = s(ia);
end

%% ===== helper: stable unique (cellstr) =====
function out = uniqueStableCell(c)
s = string(c(:));
[~, ia] = unique(s, 'stable');
out = cellstr(s(ia));
end

%% ===================== Helper: parallel pool =====================
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