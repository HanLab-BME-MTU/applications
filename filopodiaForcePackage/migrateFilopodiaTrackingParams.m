function migrateFilopodiaTrackingParams(MLPath)
%MIGRATEFILOPODIATRACKINGPARAMS  Update P3 (FilopodiaTrackingProcess) funParams
%on every movie of a MovieList from the old nearest-neighbour format
%(MaxLinkDist, LinkUseBase, ...) to the u-track format (probDim, costMatrices,
%gapCloseParam, kalmanFunctions). Run once after switching P3 to the u-track
%TrackingProcess subclass, so the stock tracking GUI opens without errors.
%
%   migrateFilopodiaTrackingParams('/path/to/movieList.mat')
%   migrateFilopodiaTrackingParams(ML)
%
% Existing tracking RESULTS are left as-is; you should re-run P3 with the
% u-track engine afterwards. Only the parameters are updated here.
% Sangyoon J. Han / 2026

if isa(MLPath,'MovieList'), ML = MLPath; else, ML = MovieList.load(MLPath,'askUser',false); end
ML.sanityCheck;

n = numel(ML.movies_); nFixed = 0;
for k = 1:n
    MD = ML.movies_{k};
    ip = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    if isempty(ip), continue; end
    proc = MD.processes_{ip};
    fp = proc.funParams_;

    % The real test is whether this process object still carries the OLD worker
    % (@trackMovieFilopodia) instead of u-track's @trackMovie. A previous
    % params-only migration may have added u-track fields to funParams while
    % leaving the old funName_/class behaviour in place, so checking funParams
    % alone is not enough.
    workerStr = '';
    try, workerStr = func2str(proc.funName_); catch, end %#ok<CTCH>
    % old worker is @trackMovieFilopodia; u-track worker is @trackMovie
    usesUtrackWorker = contains(workerStr, 'trackMovie') && ~contains(workerStr, 'trackMovieFilopodia');

    if usesUtrackWorker
        fprintf('movie %d: already u-track tracking (worker=%s), skip\n', k, workerStr);
        continue;
    end
    fprintf('movie %d: old worker=%s -> rebuilding\n', k, workerStr);

    def = FilopodiaTrackingProcess.getDefaultParams(MD);
    % preserve OutputDirectory and channel if they were set
    if isfield(fp,'OutputDirectory') && ~isempty(fp.OutputDirectory)
        def.OutputDirectory = fp.OutputDirectory;
    end

    % Replacing only funParams is not enough: a process saved by the OLD class
    % still carries the old worker handle (funName_ = @trackMovieFilopodia) and
    % the old loadChannelOutput. Rebuild the process from the current class so
    % funName_ becomes the u-track @trackMovie, then swap it into the package
    % slot in place (keeping its position and links).
    newProc = FilopodiaTrackingProcess(MD, def.OutputDirectory, def);

    % carry over outFilePaths if a result already existed
    try
        if ~isempty(proc.outFilePaths_)
            newProc.setOutFilePaths(proc.outFilePaths_);
        end
    catch
    end

    % swap old -> new via the standard API (delete then add), then rewire the
    % package slot so the package points at the rebuilt process.
    pkgIdx = MD.getPackageIndex('FilopodiaForcePackage');
    MD.deleteProcess(ip);
    MD.addProcess(newProc);
    if ~isempty(pkgIdx)
        pkg = MD.packages_{pkgIdx};
        for s = 1:numel(pkg.processes_)
            if isempty(pkg.processes_{s}) || isa(pkg.processes_{s},'FilopodiaTrackingProcess')
                pkg.setProcess(s, newProc);
            end
        end
    end
    MD.save();
    nFixed = nFixed + 1;
    fprintf('movie %d: migrated to u-track tracking (process rebuilt)\n', k);
end
try, ML.save(); catch, end %#ok<CTCH>
fprintf('Done: %d/%d movies migrated.\n', nFixed, n);
end
