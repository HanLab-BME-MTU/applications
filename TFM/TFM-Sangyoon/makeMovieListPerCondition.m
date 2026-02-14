function ML = makeMovieListPerCondition(condRoot, varargin)
% makeMovieListPerCondition
% condRoot ???? movieData.mat?? ?? MovieList(movieList.mat)? ????.
%
% ??:
%   ML = makeMovieListPerCondition("/path/to/TFMprep_bestZ/control");

ip = inputParser;
ip.addParameter('movieListName','movieList.mat', @(s)ischar(s)||isstring(s));
ip.addParameter('recursive', true, @(x)islogical(x)&&isscalar(x));
ip.parse(varargin{:});
movieListName = char(ip.Results.movieListName);
doRec = ip.Results.recursive;

condRoot = char(condRoot);
assert(exist(condRoot,'dir')==7, 'condRoot not found: %s', condRoot);

if doRec
    d = dir(fullfile(condRoot, '**', 'movieData.mat'));
else
    d = dir(fullfile(condRoot, 'movieData.mat'));
end
assert(~isempty(d), 'No movieData.mat found under: %s', condRoot);

% MovieData ??
MDs = MovieData.empty(0,1);
for k = 1:numel(d)
    mdPath = fullfile(d(k).folder, d(k).name);
    try
        md = MovieData.load(mdPath);
        MDs(end+1,1) = md; %#ok<AGROW>
    catch ME
        warning('Failed to load %s (%s). Skipping.', mdPath, ME.message);
    end
end
assert(~isempty(MDs), 'No MovieData objects could be loaded under: %s', condRoot);

% MovieList ?? (?? ??? condRoot)
ML = MovieList(MDs, condRoot);
ML.setPath(condRoot);
ML.setFilename(movieListName);
ML.sanityCheck;
ML.save;

fprintf('[makeMovieListPerCondition] Saved: %s\n', fullfile(condRoot, movieListName));
end