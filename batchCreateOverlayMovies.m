function batchCreateOverlayMovies(rootDirectory, subDirs, movieDataSubDir)

if nargin < 1 || isempty(rootDirectory)
    dataDirectory = uigetdir('', 'Select a data directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

if nargin < 2 || ~numel(subDirs) || ~iscell(subDirs)
    error('no valid subfolders provided.');
end

if nargin < 3 || isempty(movieDataSubDir)
    movieDataSubDir = subDirs{1};
end

% Get every path from rootDirectory containing the given subfolders.
nSubDirs = numel(subDirs);
paths = getDirectories(rootDirectory, nSubDirs, subDirs);

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

disp('Process all directories (Grab a coffee)...');

nMovies = numel(paths);
scrsz = get(0,'ScreenSize');
figure('Position',scrsz);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    try
        % Load movieData
        path = fullfile(paths{iMovie}, movieDataSubDir);
        filename = fullfile(path, 'movieData.mat');
        load(filename);

        % Create inputMovieInfo        
        inputMovieInfo.dir = movieData.channels(1).roiDirectory;
        filenames = dir([inputMovieInfo.dir filesep '*.tif']);
        inputMovieInfo.filename = filenames(1).name;

        % Create overlaySegment2DMovie
        if isfield(movieData,'detection') && isfield(movieData.detection, 'status') && ...
                movieData.detection.status == 1
            saveMovieInfo.dir = movieData.detection.directory;
            saveMovieInfo.filename = 'segment2DMovie.mov';
        
            load([movieData.detection.directory filesep movieData.detection.filename]);
            
            clf;
            overlaySegment2DMovie(segmentParams,[],inputMovieInfo,saveMovieInfo);
        end
        
        % Create overlayTrackedSegment2DMovie
        if isfield(movieData,'tracking') && isfield(movieData.tracking,'status') && ...
                movieData.tracking.status == 1
            saveMovieInfo.dir = movieData.tracking.directory;
            saveMovieInfo.filename = 'trackedSegment2DMovie.mov';
    
            load([movieData.tracking.directory filesep movieData.tracking.filename]);
            
            % remove tracks with lifetime == 1
            trackSEL = getTrackSEL(tracksFinal);
            tracksFinal = tracksFinal(trackSEL(:,3) >= 2);
            
            clf;
            overlayTrackedSegment2DMovie(tracksFinal,segmentParams,[],1000,inputMovieInfo,saveMovieInfo);
        end
        
    catch errMess
        disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
        disp(['Error in movie ' num2str(iMovie) ': ' errMess.message '(SKIPPING)']);
        continue;
    end
end