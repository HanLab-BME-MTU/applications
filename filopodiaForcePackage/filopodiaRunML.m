function [] = filopodiaRunML(MLPath, processesToRun, forceRun, iChanTal, useParFor)
%FILOPODIARUNML  Run the FilopodiaForcePackage on every movie in a MovieList,
%mirroring tfmRunML. Each movie is handled by filopodiaRun, which runs only
%the processes that are not already done (unless forceRun).
%
%   filopodiaRunML(MLPath)                       % run whatever is not done
%   filopodiaRunML(MLPath, 6, true)              % force-rerun P6 on all movies
%   filopodiaRunML(MLPath, 4:6, false, 1, false)
%
% Note: parfor is OFF by default because the processes write into the same
% MovieData and call MD.save; parallel saves of the same ML can collide.
% Sangyoon J. Han / 2026

if nargin < 5 || isempty(useParFor),    useParFor = false; end
if nargin < 4 || isempty(iChanTal),     iChanTal = 1; end
if nargin < 3 || isempty(forceRun),     forceRun = false; end
if nargin < 2,                          processesToRun = []; end

%% set up
if isa(MLPath,'MovieList')
    ML = MovieList.load(MLPath.getFullPath);
else
    ML = MovieList.load(MLPath,'askUser',false);
end

nMovies = numel(ML.movieDataFile_);
MDAll   = ML.movies_;

if useParFor
    parfor ii = 1:nMovies
        curMD = MDAll{ii};
        try
            filopodiaRun(curMD, processesToRun, forceRun, iChanTal);
        catch ME
            warning('Movie %d failed: %s', ii, ME.message);
        end
    end
else
    for ii = 1:nMovies
        fprintf('\n----- movie %d/%d -----\n', ii, nMovies);
        curMD = MDAll{ii};
        try
            filopodiaRun(curMD, processesToRun, forceRun, iChanTal);
        catch ME
            warning('Movie %d failed: %s', ii, ME.message);
            fprintf(2,'%s\n', getReport(ME,'extended','hyperlinks','off'));
        end
    end
end
ML.save
end
