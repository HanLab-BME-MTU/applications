%% run_TFMpipeline_allConditions_20260208.m
% Full automatic TFM pipeline:
%  0) setup TFMPackage for each MovieData in each MovieList
%  1) run TFMPackage processes: 1 (reg), 2 (displ), 4 (force)   [3 deleted]
%  2) compute strain energy by calculateMovieStrainEnergy (per MD)
%
% Requirements on path:
%   - setupTFMPackageForMovieList.m (yours)
%   - calculateMovieStrainEnergy.m (yours)
%   - MovieList/MovieData/TFMPackage classes are on path
%
% Sangyoon Han / 2026-02-08 style

clear; clc;

%% ===== USER SETTINGS =====
useParallel = true;     % set false if no PCT
maxWorkers  = 16;       % adjust to your machine
runSteps    = [1 2 4];  % processes to run in TFMPackage (3 is deleted)
doStrainEnergy = true;

% If metadata missing in MDs, fill from these defaults:
tfmParams = struct();
tfmParams.beadChan     = 2;
tfmParams.timeInterval = 120;   % sec/frame
tfmParams.pixelSize    = 325;   % nm/pixel (update if needed)
tfmParams.YoungModulus = 5000;  % Pa
tfmParams.regParam     = 1e-5;

% Displacement params (optional)
tfmParams.alpha = 0.05;
tfmParams.minCorLength = 15;
tfmParams.addNonLocMaxBeads = true;
tfmParams.maxFlowSpeed = 20;
tfmParams.highRes = true;
tfmParams.useGrid = true;
tfmParams.mode = 'accurate';
tfmParams.noFlowOutwardOnBorder = 1;
tfmParams.trackSuccessively = false;

%% ===== MovieList paths (EDIT THESE) =====
MLpaths = {
  "/mnt/nas/shear_stress/Inflammation/02_8_2026_control_1%DMSO/02_8_2026_control_1%DMSO_1/TFMprep_bestZ/control/MovieList/movieList.mat"
  "/mnt/nas/shear_stress/Inflammation/02_8_2026_TNFa10ng/02_8_2026_TNFa10ng_noDMSO_1/TFMprep_bestZ/TNFa/MovieList/movieList.mat"
  "/mnt/nas/shear_stress/Inflammation/02_8_2026_BB100um+TNF10ngml/02_8_2026_BB100um+TNF10ngml_1/TFMprep_bestZ/TNFa_BBS/MovieList/movieList.mat"
  "/mnt/nas/shear_stress/Inflammation/02_8_2026_BB100um/02_8_2026_BB100um_1/TFMprep_bestZ/BBS/MovieList/movieList.mat"
};

%% ===== MAIN =====
if useParallel
    try
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local', maxWorkers);
        end
    catch ME
        warning('Parallel pool failed (%s). Falling back to serial.', ME.message);
        useParallel = false;
    end
end

for ci = 1:numel(MLpaths)
    mlPath = char(MLpaths{ci});
    fprintf('\n====================================================\n');
    fprintf('Condition %d/%d\nML: %s\n', ci, numel(MLpaths), mlPath);
    fprintf('====================================================\n');

    assert(exist(mlPath,'file')==2, 'movieList not found: %s', mlPath);

    % 0) Setup TFMPackage parameters
    fprintf('[0] Setting up TFMPackage...\n');
    setupTFMPackageForMovieList(mlPath, tfmParams);

    % Reload MovieList after setup
    ML = MovieList.load(mlPath, 'askUser', false);
    nMovies = numel(ML.movieDataFile_);
    fprintf('[Info] nMovies=%d\n', nMovies);

    % Collect MD objects
    MDs = cell(nMovies,1);
    for i=1:nMovies
        MDs{i} = ML.getMovie(i);
    end

    % 1) Run TFMPackage steps
    fprintf('[1] Running TFMPackage steps %s ...\n', mat2str(runSteps));

    tfmRunML(ML,useParallel);

    fprintf('[Done] %s\n', mlPath);
end

disp('All conditions finished.');

%% ===== Local function =====
function runTFMStepsForOneMD(MD, steps)
% Robust runner for TFMPackage processes.
% Assumes setupTFMPackageForMovieList already configured parameters.

try
    iPack = MD.getPackageIndex('TFMPackage');
    assert(~isempty(iPack), 'No TFMPackage in MD: %s', MD.outputDirectory_);
    pack = MD.getPackage(iPack);

    for s = steps
        if numel(pack.processes_) < s || isempty(pack.processes_{s})
            error('TFMPackage missing process %d in %s', s, MD.outputDirectory_);
        end
        proc = pack.processes_{s};

        % Some Process classes expose run(); others might use runFunction().
        if ismethod(proc,'run')
            proc.run();
        elseif ismethod(proc,'runFunction')
            proc.runFunction();
        else
            error('Process %d has no run method (%s)', s, class(proc));
        end
    end

catch ME
    warning('TFM steps failed: %s', ME.message);
end
end