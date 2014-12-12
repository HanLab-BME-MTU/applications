function batchCreateVinculinMovies(rootDirectory)

if nargin < 1 || isempty(rootDirectory)
    dataDirectory = uigetdir('', 'Select a data directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

% Get every path from rootDirectory containing ch488 and ch560 folders
paths = getDirectories(rootDirectory, 2, {'ch488','ch560'});

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

movieDataSubDir = ['ch488' filesep 'analysis'];

nMovies = numel(paths);
scrsz = get(0,'ScreenSize');
figure('Position',scrsz);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    try
        %% ------- Load movieData ----------- %%
        path = fullfile(paths{iMovie}, movieDataSubDir);
        filename = fullfile(path, 'movieData.mat');
        load(filename);

        %% ------- Create inputMovieInfo ----------- %%
        inputMovieInfo.dir = fullfile(movieData.imageDirectory, movieData.channelDirectory{1});
        filenames = dir([inputMovieInfo.dir filesep '*.tif']);
        inputMovieInfo.filename = filenames(1).name;

        %% ------- Create segment detection movie -------- %%

        if checkMovieSegmentDetection(movieData)
            saveMovieInfo.dir = movieData.segmentDetection.directory;
            saveMovieInfo.filename = 'segment2DMovie.mov';
        
            load([movieData.segmentDetection.directory filesep movieData.segmentDetection.filename]);
            
            clf;
            overlaySegment2DMovie(segmentParams,[],inputMovieInfo,saveMovieInfo);
        end
        
        %% ------- Create segment tracking movie --------%%
        
        if checkMovieSegmentTracking(movieData)
            saveMovieInfo.dir = movieData.segmentTracking.directory;
            saveMovieInfo.filename = 'trackedSegment2DMovie.mov';
    
            load([movieData.segmentTracking.directory filesep movieData.segmentTracking.filename]);
            
            % remove tracks with lifetime == 1
            %trackSEL = getTrackSEL(tracksFinal);
            %tracksFinal = tracksFinal(trackSEL(:,3) >= 2);
            
            clf;
            overlayTrackedSegment2DMovie(tracksFinal,segmentParams,[],1000,inputMovieInfo,saveMovieInfo);
        end
        
    catch errMess
        disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
        disp(['Error in movie ' num2str(iMovie) ': ' errMess.message '(SKIPPING)']);
        continue;
    end
end