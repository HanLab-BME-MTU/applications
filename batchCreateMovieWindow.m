function batchCreateMovieWindow(varargin)

if nargin == 1
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

% Get all subdirectories containing Actin & TM.
paths = getDirectories(rootDirectory);
nMovies = numel(paths);

movieData(1:nMovies) = struct('analysisDirectory', '',...
    'imageDirectory', [],...
    'nImages', [],...
    'channelDirectory', [],...
    'pixelSize_nm', [],...
    'timeInterval_s', []);

% Set up parameters
dContour = 15; % ~ 1um
dWin = 10;
iStart = 1;
iEnd = 2;
winMethod = 'e';

for iMovie = 1:numel(movieData)
    
    % STEP 1: Create the initial movie data
    movieData(iMovie).analysisDirectory = [paths{iMovie} filesep 'windowAnalysis'];
    
    % Be sure to handle 'Actin' lower case properly
    movieData(iMovie).imageDirectory = [paths{iMovie} filesep 'Actin' filesep 'crop'];
    if ~exist(movieData(iMovie).imageDirectory, 'dir')
        disp(['Movie ' num2str(iMovie) '/' num2str(nMovies) ': unable to find image directory.']);
        continue;
    end
    
    movieData(iMovie).nImages = numel(dir([movieData(iMovie).imageDirectory filesep '*.tif']));
    movieData(iMovie).channelDirectory = {''};
    % TODO: read the pixel size from FSM center parameter file.
    movieData(iMovie).pixelSize_nm = 67;
    % TODO: read the time interval from FSM center parameter file.
    movieData(iMovie).timeInterval_s = 10;
    
    % Get the mask directory
    movieData(iMovie).masks.directory = [paths{iMovie} filesep 'Actin' filesep 'cell_mask'];
    if ~exist(movieData(iMovie).masks.directory, 'dir')
        disp(['Movie ' num2str(iMovie) '/' num2str(nMovies) ': unable to find mask directory.']);
        continue;
    end
    movieData(iMovie).masks.n = numel(dir([movieData(iMovie).masks.directory filesep '*.tif']));
    
    % STEP 2: Get the contour
    try
        movieData(iMovie) = getMovieContours(movieData(iMovie), 0:dContour:500, 0, 0, ['contours_'  num2str(dContour) 'pix.mat']);
    catch errMess
        disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
        currMovie.contours.error = errMess;
        currMovie.contours.status = 0;
    end
    
    % STEP 3: Calculate protusion
    try
        %Indicate that protrusion was started
        movieData(iMovie).protrusion.status = 0;

        %Set up inputs for sams protrusion calc

        %Check and convert OS for linux
        movieData(iMovie) = setupMovieData(movieData(iMovie));

        handles.batch_processing = 1;
        handles.directory_name = [movieData(iMovie).masks.directory];
        handles.result_directory_name = [movieData(iMovie).masks.directory];
        handles.FileType = '*.tif';
        handles.timevalue = movieData(iMovie).timeInterval_s;
        handles.resolutionvalue = movieData(iMovie).pixelSize_nm;
        handles.segvalue = 30;
        handles.dl_rate = 50;
        
        %run it
        [OK,handles] = protrusionAnalysis(handles);

        if ~OK
            movieData(iMovie).protrusion.status = 0;
        else
            if isfield(currMovie.protrusion,'error')
                movieData(iMovie).protrusion = rmfield(movieData(iMovie).protrusion,'error');
            end
            movieData(iMovie).protrusion.directory = [movieData(iMovie).masks.directory];
            movieData(iMovie).protrusion.fileName = 'protrusion.mat';
            movieData(iMovie).protrusion.nfileName = 'normal_matrix.mat';
            movieData(iMovie).protrusion.status = 1;
        end
    catch errMess
        disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
        movieData(iMovie).protrusion.error = errMess;
        movieData(iMovie).protrusion.status = 0;
    end
    
    % STEP 4: Create windows    
    windowString = [num2str(dContour) 'by' num2str(dWin) 'pix_' num2str(iStart) '_' num2str(iEnd)];
     
    try
        movieData(iMovie) = setupMovieData(movieData(iMovie));
            
        disp(['Windowing movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
        movieData(iMovie) = getMovieWindows(movieData(iMovie),winMethod,dWin,[],iStart,iEnd,[],[],['windows_' winMethod '_' windowString '.mat']);
        
        if isfield(currMovie.windows,'error')
            movieData(iMovie).windows = rmfield(movieData(iMovie).windows,'error');
        end
    
    catch errMess
        disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
        movieData(iMovie).windows.error = errMess;
        movieData(iMovie).windows.status = 0;
    end
end
end