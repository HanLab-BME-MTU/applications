function batchCreateMovieWindow(varargin)

if nargin >= 1
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

forceRun = 0;
if nargin == 2
    forceRun = varargin{2};
end

% Get all subdirectories containing Actin & TM.
paths = getDirectories(rootDirectory);
nMovies = numel(paths);

movieData = cell(nMovies, 1);

% Set up parameters
dContour = 15; % ~ 1um
dWin = 10;
iStart = 1;
iEnd = 10;
winMethod = 'e';

for iMovie = 1:numel(movieData)

    %Get current movie data for readability
    currMovie = movieData{iMovie};
    
    % STEP 1: Create the initial movie data
    currMovie.analysisDirectory = [paths{iMovie} filesep 'windowAnalysis'];
    
    % Be sure to handle 'Actin' lower case properly
    currMovie.imageDirectory = [paths{iMovie} filesep 'Actin' filesep 'crop'];
    if ~exist(currMovie.imageDirectory, 'dir')
        disp(['Movie ' num2str(iMovie) '/' num2str(nMovies) ': unable to find image directory.']);
        continue;
    end

    currMovie.nImages = numel(dir([currMovie.imageDirectory filesep '*.tif']));
    currMovie.channelDirectory = {''};
    % TODO: read the pixel size from FSM center parameter file.
    currMovie.pixelSize_nm = 67;
    % TODO: read the time interval from FSM center parameter file.
    currMovie.timeInterval_s = 10;

    % Get the mask directory
    currMovie.masks.directory = [paths{iMovie} filesep 'Actin' filesep 'edge' filesep 'cell_mask'];
    if ~exist(currMovie.masks.directory, 'dir')
        disp(['Movie ' num2str(iMovie) '/' num2str(nMovies) ': unable to find mask directory.']);
        continue;
    end
    currMovie.masks.n = numel(dir([currMovie.masks.directory filesep '*.tif']));

    % Update from already saved movieData
    if exist([currMovie.analysisDirectory filesep 'movieData.mat'], 'file') && ~forceRun
        currMovie = refreshMovieData(currMovie);
    end

    % STEP 2: Get the contour
    if ~isfield(currMovie,'contours') || ~isfield(currMovie.contours,'status') || currMovie.contours.status ~= 1 || forceRun
        try
            currMovie = getMovieContours(currMovie, 0:dContour:500, 0, 0, ['contours_'  num2str(dContour) 'pix.mat']);
        catch errMess
            disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
            currMovie.contours.error = errMess;
            currMovie.contours.status = 0;
        end
    end

    % STEP 3: Calculate protusion
    if ~isfield(currMovie,'protrusion') || ~isfield(currMovie.protrusion,'status') || currMovie.protrusion.status ~= 1 || forceRun
        try
            %Indicate that protrusion was started
            currMovie.protrusion.status = 0;

            %Set up inputs for sams protrusion calc

            %Check and convert OS for linux
            currMovie = setupMovieData(currMovie);

            handles.batch_processing = 1;
            handles.directory_name = [currMovie.masks.directory];
            handles.result_directory_name = [currMovie.masks.directory];
            handles.FileType = '*.tif';
            handles.timevalue = currMovie.timeInterval_s;
            handles.resolutionvalue = currMovie.pixelSize_nm;
            handles.segvalue = 30;
            handles.dl_rate = 50;

            %run it
            [OK,handles] = protrusionAnalysis(handles);

            if ~OK
                currMovie.protrusion.status = 0;
            else
                if isfield(currMovie.protrusion,'error')
                    currMovie.protrusion = rmfield(currMovie.protrusion,'error');
                end
                currMovie.protrusion.directory = [currMovie.masks.directory];
                currMovie.protrusion.fileName = 'protrusion.mat';
                currMovie.protrusion.nfileName = 'normal_matrix.mat';
                currMovie.protrusion.status = 1;
            end
            
            updateMovieData(currMovie);
            
        catch errMess
            disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
            currMovie.protrusion.error = errMess;
            currMovie.protrusion.status = 0;
        end
    end

    % STEP 4: Create windows
    windowString = [num2str(dContour) 'by' num2str(dWin) 'pix_' num2str(iStart) '_' num2str(iEnd)];

    if ~isfield(currMovie,'windows') || ~isfield(currMovie.windows,'status')  || currMovie.windows.status ~= 1 || forceRun
        try
            currMovie = setupMovieData(currMovie);

            disp(['Windowing movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieWindows(currMovie,winMethod,dWin,[],iStart,iEnd,[],[],['windows_' winMethod '_' windowString '.mat']);

            if isfield(currMovie.windows,'error')
                currMovie.windows = rmfield(currMovie.windows,'error');
            end

        catch errMess
            disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
            currMovie.windows.error = errMess;
            currMovie.windows.status = 0;
        end
    end

    % STEP 5: Create the window labels;
    if ~isfield(currMovie, 'labels') || ~isfield(currMovie.labels, 'status') || currMovie.labels.status ~= 1 || forceRun
        try
            currMovie = setupMovieData(currMovie);

            disp(['Labeling movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            
            currMovie = getMovieLabels(currMovie);

            if isfield(currMovie.labels,'error')
                currMovie.labels = rmfield(currMovie.labels,'error');
            end

        catch errMess
            disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
            currMovie.labels.error = errMess;
            currMovie.labels.status = 0;
        end            
    end

    % STEP 6: Split the windows into different files.
    if ~isfield(currMovie, 'windows') || ~isfield(currMovie.windows, 'splitted') || currMovie.windows.splitted ~= 1
        splitWindowFrames(currMovie, [currMovie.analysisDirectory filesep 'windows']);
        currMovie.windows.splitted = 1;
    end

    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMov) ': ' errMess.message],mfileName());
    end
    
    movieData{iMovie} = currMovie;
end

end