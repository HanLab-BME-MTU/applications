%% runFilopodiaForceBatch.m
% Load the movieList.mat in each condition folder and run the full
% FilopodiaForcePackage (P1-P6) on every movie, using the default settings
% from runFilopodiaForce_debug.
%
% Edit rootDir and condFolders to match your layout, then run.
% After this finishes, use filopodiaForceBatch.m to compare conditions.
% Sangyoon J. Han / 2026

%% ===================== CONFIG =====================
rootDir = '/mnt/nas/Collaborations/Celine_San_Diego/TFM_Comparison_Nov-Apr2026';

% condition folders to process (each must contain movieList.mat)
condFolders = { ...
    'WT-Soft', ...
    'KO-Soft', ...
    '1YA-Soft', ...
    'WT-Stiff', ...
    'KO-Stiff' };
% NOTE: add/remove to match what you want to run (e.g. include 1YA-Stiff if present)

iChanTal  = 1;      % talin-GFP channel
overwrite = false;  % true = re-run processes even if outputs already exist
mlName    = 'movieList.mat';

%% ===================== RUN =====================
nCond = numel(condFolders);
for c = 1:nCond
    condDir = fullfile(rootDir, condFolders{c});
    mlPath  = fullfile(condDir, mlName);
    fprintf('\n========== Condition %d/%d: %s ==========\n', c, nCond, condFolders{c});

    if exist(mlPath,'file')~=2
        warning('No %s in %s, skipping condition.', mlName, condDir);
        continue;
    end

    ML = MovieList.load(mlPath);
    ML.sanityCheck;
    nMov = numel(ML.movies_);
    fprintf('Loaded MovieList with %d movies.\n', nMov);

    for k = 1:nMov
        fprintf('\n----- [%s] movie %d/%d -----\n', condFolders{c}, k, nMov);
        MD = ML.movies_{k};
        try
            runFilopodiaForceOnMD(MD, iChanTal, 'Overwrite', overwrite);
            MD.save();
        catch ME
            warning('Movie %d in %s failed: %s', k, condFolders{c}, ME.message);
            fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks','off'));
        end
    end

    % persist process registration on the MovieList movies
    try, ML.save(); catch, end %#ok<NOSEMI,CTCH>
end

fprintf('\n========== runFilopodiaForceBatch done ==========\n');
fprintf('Next: run filopodiaForceBatch.m (choose the condition folders) to compare.\n');
