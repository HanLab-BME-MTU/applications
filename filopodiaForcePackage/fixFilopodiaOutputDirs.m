function fixFilopodiaOutputDirs(ML)
%FIXFILOPODIAOUTPUTDIRS  Repair OutputDirectory for P1-P6 Filopodia processes
%across all movies in a MovieList. Also removes duplicate process instances
%that can accumulate from repeated migrations.
%
%   fixFilopodiaOutputDirs(ML)
%   fixFilopodiaOutputDirs('path/to/movieList.mat')
%
% The "apply to all movies" operation in the settings GUI copies funParams
% including OutputDirectory from the source movie to all others, causing
% every movie to write results to the same folder. This function rebuilds
% each process's OutputDirectory from the movie's own outputDirectory_.
% Sangyoon J. Han / 2026

if ischar(ML), ML = MovieList.load(ML,'askUser',false); end
ML.sanityCheck;

procClasses = {
    'FilopodiaSegmentationProcess', 'FilopodiaSegmentation';
    'FilopodiaDetectionProcess',    'FilopodiaDetection';
    'FilopodiaClassificationProcess','FilopodiaClassification';
    'FilopodiaSamplingProcess',     'FilopodiaSampling';
    'FilopodiaStatisticsProcess',   'FilopodiaStatistics';
};
pkgSubdir = 'FilopodiaForcePackage';

for k = 1:numel(ML.movies_)
    MD = ML.movies_{k};
    mdRoot = MD.outputDirectory_;
    changed = false;

    for c = 1:size(procClasses,1)
        cls   = procClasses{c,1};
        subdir = procClasses{c,2};
        correctDir = fullfile(mdRoot, pkgSubdir, subdir);

        % find ALL process indices of this class (may be duplicates)
        idxs = [];
        for i = 1:numel(MD.processes_)
            p = MD.processes_{i};
            if ~isempty(p) && isa(p, cls)
                idxs(end+1) = i; %#ok<AGROW>
            end
        end
        if isempty(idxs), continue; end

        % remove duplicates: keep the last one (most recent run) and delete rest
        if numel(idxs) > 1
            fprintf('movie %d: removing %d duplicate %s\n', k, numel(idxs)-1, cls);
            % keep last index, delete others
            for di = idxs(1:end-1)
                try, MD.deleteProcess(MD.processes_{di}); catch, end
            end
            % re-find after deletion
            idxs = [];
            for i = 1:numel(MD.processes_)
                p = MD.processes_{i};
                if ~isempty(p) && isa(p, cls)
                    idxs(end+1) = i; %#ok<AGROW>
                end
            end
        end

        if isempty(idxs), continue; end
        proc = MD.processes_{idxs(1)};

        % fix OutputDirectory if it doesn't belong to this movie
        if ~isfield(proc.funParams_,'OutputDirectory') || ...
           ~startsWith(proc.funParams_.OutputDirectory, mdRoot)
            fp = proc.funParams_;
            fp.OutputDirectory = correctDir;
            try, proc.setPara(fp); catch, proc.funParams_ = fp; end
            fprintf('movie %d: fixed %s -> %s\n', k, cls, correctDir);
            changed = true;
        end
    end

    if changed
        MD.save();
        fprintf('movie %d: saved\n', k);
    else
        fprintf('movie %d: OK (no changes needed)\n', k);
    end
end
try, ML.save(); catch, end
fprintf('Done.\n');
end
