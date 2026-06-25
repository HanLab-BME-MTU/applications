function removeDuplicateFilopodiaProcesses(ML)
%REMOVEDUPLICATEFILOPODIAPROCESSES  Remove duplicate Filopodia process
%instances from every movie in a MovieList. Duplicates accumulate when
%settings are applied across movies or migration scripts run multiple times.
%Keeps the LAST instance of each class (most recent run) and removes others.
%Also fixes OutputDirectory to each movie's own path.
%
%   removeDuplicateFilopodiaProcesses(ML)
%   removeDuplicateFilopodiaProcesses('path/to/movieList.mat')
%
% Sangyoon J. Han / 2026

if ischar(ML), ML = MovieList.load(ML,'askUser',false); end
ML.sanityCheck;

filoClasses = {
    'FilopodiaSegmentationProcess'
    'FilopodiaDetectionProcess'
    'FilopodiaClassificationProcess'
    'FilopodiaSamplingProcess'
    'FilopodiaStatisticsProcess'
    'FilopodiaTrackingProcess'
};

pkgSubdir = 'FilopodiaForcePackage';
procSubdirs = containers.Map(filoClasses, {
    'FilopodiaSegmentation'
    'FilopodiaDetection'
    'FilopodiaClassification'
    'FilopodiaSampling'
    'FilopodiaStatistics'
    'FilopodiaTracking'
});

for k = 1:numel(ML.movies_)
    MD = ML.movies_{k};
    mdRoot = MD.outputDirectory_;
    changed = false;

    for c = 1:numel(filoClasses)
        cls = filoClasses{c};

        % find ALL indices of this class
        idxs = [];
        for i = 1:numel(MD.processes_)
            p = MD.processes_{i};
            if ~isempty(p) && strcmp(class(p), cls)
                idxs(end+1) = i; %#ok<AGROW>
            end
        end
        if isempty(idxs), continue; end

        if numel(idxs) > 1
            fprintf('movie %d: %d duplicates of %s -> keeping index %d, removing %s\n', ...
                k, numel(idxs)-1, cls, idxs(end), mat2str(idxs(1:end-1)));
            % keep the last one; delete earlier duplicates in reverse order
            % so each deletion doesn't shift the indices of remaining targets
            toRemove = sort(idxs(1:end-1), 'descend');
            for di = toRemove
                try
                    proc_del = MD.processes_{di};
                    MD.deleteProcess(proc_del);
                catch
                    % deleteProcess may fail on invalid handles; try direct removal
                    try
                        if di <= numel(MD.processes_)
                            MD.processes_(di) = [];
                        end
                    catch
                    end
                end
            end
            changed = true;
        end

        % re-find the surviving process after any deletions
        idxs = [];
        for i = 1:numel(MD.processes_)
            p = MD.processes_{i};
            if ~isempty(p) && strcmp(class(p), cls)
                idxs(end+1) = i; %#ok<AGROW>
            end
        end

        if isempty(idxs), continue; end
        proc = MD.processes_{idxs(1)};
        if isempty(proc) || ~isvalid(proc), continue; end

        % fix OutputDirectory if it doesn't belong to this movie
        if isKey(procSubdirs, cls)
            correctDir = fullfile(mdRoot, pkgSubdir, procSubdirs(cls));
            fp = proc.funParams_;
            if ~isfield(fp,'OutputDirectory') || ...
               ~startsWith(fp.OutputDirectory, mdRoot)
                fp.OutputDirectory = correctDir;
                try, proc.setPara(fp); catch, proc.funParams_ = fp; end
                fprintf('movie %d: fixed OutputDirectory for %s\n', k, cls);
                changed = true;
            end
        end
    end

    % rewire package slots after process list changed
    try
        pkgIdx = MD.getPackageIndex('FilopodiaForcePackage');
        if ~isempty(pkgIdx)
            pkg = MD.packages_{pkgIdx};
            pkg.sanityCheck();
        end
    catch
    end

    if changed
        MD.save();
        fprintf('movie %d: saved\n', k);
    else
        fprintf('movie %d: OK\n', k);
    end
end
try, ML.save(); catch, end
fprintf('Done.\n');
end
