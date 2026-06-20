%% rerunFilopodiaP6.m
% Re-run ONLY Process 6 (FilopodiaStatisticsProcess) on movies where it is
% missing or did not finish successfully. P1-P5 are left untouched.
%
% Use this after fixing computeMovieFilopodiaStats.m, to recover movies that
% threw in P6 without re-running the whole pipeline.
%
% IMPORTANT: make sure the corrected computeMovieFilopodiaStats.m is on the
% MATLAB path (run "rehash path" or "clear functions" first if you just
% replaced the file).
% Sangyoon J. Han / 2026

%% ===================== CONFIG =====================
rootDir = '/mnt/nas/Collaborations/Celine_San_Diego/TFM_Comparison_Nov-Apr2026';
condFolders = { ...
    'WT-Soft', ...
    'KO-Soft', ...
    '1YA-Soft', ...
    'WT-Stiff', ...
    'KO-Stiff' };
iChanTal = 1;
mlName   = 'movieList.mat';
forceAll = false;   % true = rerun P6 on every movie regardless of state

%% ===================== RUN =====================
clear functions   %#ok<CLFUNC>  % ensure fixed worker is picked up
nCond = numel(condFolders);
nFixed = 0; nSkipped = 0; nFailed = 0;

for c = 1:nCond
    mlPath = fullfile(rootDir, condFolders{c}, mlName);
    if exist(mlPath,'file')~=2
        warning('No %s in %s, skipping.', mlName, condFolders{c}); continue;
    end
    ML = MovieList.load(mlPath); ML.sanityCheck;
    nMov = numel(ML.movies_);
    fprintf('\n===== %s (%d movies) =====\n', condFolders{c}, nMov);

    for k = 1:nMov
        MD = ML.movies_{k};
        iP6 = MD.getProcessIndex('FilopodiaStatisticsProcess',1,0);

        % does P6 already exist and look complete?
        done = false;
        if ~isempty(iP6)
            statProc = MD.processes_{iP6};
            okSucc = false;
            try, okSucc = ~isempty(statProc.success_) && statProc.success_; catch, end %#ok<CTCH>
            okFile = false;
            try
                f = statProc.outFilePaths_{1,iChanTal};
                if isempty(f)||exist(f,'file')~=2
                    f = fullfile(statProc.funParams_.OutputDirectory,'filoStats.mat');
                end
                okFile = exist(f,'file')==2;
            catch
            end
            done = okSucc && okFile;
        end

        if done && ~forceAll
            fprintf('  movie %d: P6 OK, skip\n', k); nSkipped = nSkipped + 1; continue;
        end

        % make sure P5 (its input) is present before attempting P6
        iP5 = MD.getProcessIndex('FilopodiaSamplingProcess',1,0);
        p5ok = false;
        if ~isempty(iP5)
            sp = MD.processes_{iP5};
            try
                f5 = sp.outFilePaths_{1,iChanTal};
                if isempty(f5)||exist(f5,'file')~=2
                    f5 = fullfile(sp.funParams_.OutputDirectory,'filoSamples.mat');
                end
                p5ok = exist(f5,'file')==2;
            catch
            end
        end
        if ~p5ok
            warning('  movie %d: P5 output missing, cannot run P6. Run P1-P5 first.', k);
            nFailed = nFailed + 1; continue;
        end

        % add P6 if not registered
        if isempty(iP6)
            MD.addProcess(FilopodiaStatisticsProcess(MD));
            iP6 = MD.getProcessIndex('FilopodiaStatisticsProcess',1,0);
            statProc = MD.processes_{iP6};
            pp6 = statProc.funParams_;
            pp6.ChannelIndex        = iChanTal;
            pp6.MinLifetimeForStats = 3;
            pp6.MakeFigures         = true;
            pp6.OutputDirectory     = fullfile(MD.outputDirectory_, ...
                'FilopodiaForcePackage','FilopodiaStatistics');
            statProc.setPara(pp6);
        end

        fprintf('  movie %d: running P6...\n', k);
        try
            statProc.run();
            MD.save();
            nFixed = nFixed + 1;
        catch ME
            warning('  movie %d: P6 failed again: %s', k, ME.message);
            fprintf(2,'%s\n', getReport(ME,'extended','hyperlinks','off'));
            nFailed = nFailed + 1;
        end
    end
    try, ML.save(); catch, end %#ok<NOSEMI,CTCH>
end

fprintf('\n===== rerunFilopodiaP6 done: %d fixed, %d already OK, %d failed =====\n', ...
    nFixed, nSkipped, nFailed);
