function movieData = getMoviePairTracks(movieData,maxDistance,...
    minOverlap,bandWidth,minDistance,alpha,batchMode)

% maxDistance:     maximum distance between track pair candidates (t1,t2) such
%                  that max(dist(t1,t2)) < maxDistance.
%
% minOverlap:      is the minimal number of frame 2 tracks need to live together
%                  to be considered as potential pair of tracks. Default is 1.
%
% bandWidth:       distance in nanometers away from the cell edge where pair of
%                  feature needs to be evaluate whether it follows the Actin
%                  flow.
%
% minDistance:     minimum distance between 2 track pair candidates in the
%                  first 'bandWidth' nanometers away from cell edge.
%
% alpha:           quantile of PDF tail. Default is 0.05.

% Parse input parameters

if nargin < 3 || isempty(minOverlap)
    minOverlap = 1;
end

if nargin < 4 || isempty(bandWidth)
    bandWidth = 1000;
end

if nargin < 5 || isempty(alpha)
    alpha = 0.05;
end

assert(checkMovieDistanceTransform(movieData));
assert(checkMovieParticleDetection(movieData));
assert(checkMovieParticleTracking(movieData));

% BEGIN
movieData.pairTracks.status = 0;

% Create output directory
movieData.pairTracks.directory = fullfile(movieData.analysisDirectory, 'pairTracks');

if ~exist(movieData.pairTracks.directory, 'dir')
    mkdir(movieData.pairTracks.directory);
end

% Get all track parameters
[allTrackParams, tFirst, lifetime] = getAllTrackParams(movieData);

% Get initial pair track candidates
[E, tOverlapFirst, overlap, pFirst1, pFirst2] = ...
    getInitialPairTracks(movieData, allTrackParams, tFirst, lifetime, ...
    maxDistance, minOverlap, bandWidth, minDistance);

% END
movieData.pairTracks.status = 1;
