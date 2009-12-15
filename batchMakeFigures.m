function batchMakeFigures(varargin)

if nargin >= 1 && ~isempty(varargin{1})
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

if nargin >= 2 && ~isempty(varargin{2})
    subDirNames = varargin{2};
else
    error('No subdir specified.');
end

if nargin >= 3 && ~isempty(varargin{3})
    forceRun = varargin{3};
else
    forceRun = zeros(7, 1);
end

if nargin >= 4 && ~isempty(varargin{4})
    batchMode = varargin{4};
else
    batchMode = 1;
end

% Get every path from rootDirectory containing actin subfolder.
paths = getDirectories(rootDirectory, 2, subDirNames, ...
    @(x) exist([x filesep 'lastProjSettings.mat'], 'file'));

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

disp('Process all directories (Grab a coffee)...');

nMovies = numel(paths);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    path = paths{iMovie};
   
    currMovie = movieData{iMovie};
    
    % STEP 1: Create movieData
    currMovie.analysisDirectory = [paths{iMovie} filesep 'windowAnalysis'];
    content = dir(path);
    content = {content.name};
    ind = cellfun(@(x) any(strcmpi(x, subDirNames)), content);
    if ~nnz(ind)
        disp([movieName ': Unable to find the FSM directory. (SKIPPING)']);
        continue;
    end    
    %fsmFolderName = content{find(ind, 1, 'first')};
    currMovie.fsmDirectory = cellfun(@(x) [path filesep x], content(ind), ...
        'UniformOutput', false);
    % Add these 2 fields to be compliant with Hunter's check routines:
    currMovie.imageDirectory = [currMovie.fsmDirectory{1} filesep 'crop'];
    currMovie.channelDirectory = {''};
    currMovie.nImages = numel(dir([currMovie.imageDirectory filesep '*.tif']));
    %load([path filesep fsmFolderName filesep 'fsmPhysiParam.mat']);
    load([currMovie.fsmDirectory{1} filesep 'fsmPhysiParam.mat']);
    currMovie.pixelSize_nm = fsmPhysiParam.pixelSize;
    currMovie.timeInterval_s = fsmPhysiParam.frameInterval;
    clear fsmPhysiParam;

    currMovie.masks.directory = [currMovie.fsmDirectory{1} filesep 'edge' filesep 'cell_mask'];
    if ~exist(currMovie.masks.directory, 'dir')
        disp([movieName ': unable to find mask directory. (SKIPPING)']);
        continue;
    end
    currMovie.masks.channelDirectory = {''};
    currMovie.masks.n = numel(dir([currMovie.masks.directory filesep '*.tif']));
    currMovie.masks.status = 1;
    
    if exist([currMovie.analysisDirectory filesep 'movieData.mat'], 'file') && ~forceRun(1)
       currMovie = load([currMovie.analysisDirectory filesep 'movieData.mat']);
       currMovie = currMovie.movieData;
    end
    % Temporary update of existing movieData:
%     if exist([currMovie.analysisDirectory filesep 'movieData.mat'], 'file') && forceRun(1)
%         currMovieOld = load([currMovie.analysisDirectory filesep 'movieData.mat']);
%         currMovieOld = currMovieOld.movieData;
%         % Update new fields.
%         currMovieOld.fsmDirectory = currMovie.fsmDirectory;
%         currMovieOld.imageDirectory = currMovie.imageDirectory;
%         currMovieOld.masks.directory = currMovie.masks.directory;
%         % Save it
%         updateMovieData(currMovieOld);
%         currMovie = currMovieOld;
%     end

    % STEP 2: Get contours
    dContour = 1000 / currMovie.pixelSize_nm; % ~ 1um
    dWin = 2000 / currMovie.pixelSize_nm; % ~ 2um
    iStart = 2;
    iEnd = 4;
    winMethod = 'c';    
    
    if ~isfield(currMovie,'contours') || ~isfield(currMovie.contours,'status') || ...
            currMovie.contours.status ~= 1 || forceRun(2)
        try
            disp(['Get contours of movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieContours(currMovie, 0:dContour:500, 0, 1, ...
                ['contours_'  num2str(dContour) 'pix.mat'], batchMode);
        catch errMess
            disp(['Error in movie ' num2str(iMovie) ': ' errMess.message]);
            currMovie.contours.error = errMess;
            currMovie.contours.status = 0;
            continue;
        end
    end
    
    % STEP 3: Calculate protrusion
    if ~isfield(currMovie,'protrusion') || ~isfield(currMovie.protrusion,'status') || ...
            currMovie.protrusion.status ~= 1 || forceRun(3)
        try
            currMovie.protrusion.status = 0;

            currMovie = setupMovieData(currMovie);

            handles.batch_processing = batchMode;
            handles.directory_name = [currMovie.masks.directory];
            handles.result_directory_name = [currMovie.masks.directory];
            handles.FileType = '*.tif';
            handles.timevalue = currMovie.timeInterval_s;
            handles.resolutionvalue = currMovie.pixelSize_nm;
            handles.segvalue = 30;
            handles.dl_rate = 30;

            %run it
            [OK,handles] = protrusionAnalysis(handles);

            if ~OK
                currMovie.protrusion.status = 0;
            else
                if isfield(currMovie.protrusion,'error')
                    currMovie.protrusion = rmfield(currMovie.protrusion,'error');
                end
                
                currMovie.protrusion.directory = [currMovie.masks.directory filesep ...
                    'analysis_dl' num2str(handles.dl_rate)];
                
                currMovie.protrusion.fileName = 'protrusion.mat';
                currMovie.protrusion.nfileName = 'normal_matrix.mat';
                currMovie.protrusion.status = 1;
            end
            
            updateMovieData(currMovie);
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.protrusion.error = errMess;
            currMovie.protrusion.status = 0;
            continue;
        end
    end
    
    % STEP 4: Create windowing
    windowString = [num2str(dContour) 'by' num2str(dWin) 'pix_' ...
                num2str(iStart) '_' num2str(iEnd)];
            
    if ~isfield(currMovie,'windows') || ~isfield(currMovie.windows,'status')  || ...
            currMovie.windows.status ~= 1 || forceRun(4)
        try
            currMovie = setupMovieData(currMovie);

            disp(['Windowing movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieWindows(currMovie,winMethod,dWin,[],iStart,iEnd,[],[],...
                ['windows_' winMethod '_' windowString '.mat'], batchMode);
            
            if isfield(currMovie.windows,'error')
                currMovie.windows = rmfield(currMovie.windows,'error');
            end
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.windows.error = errMess;
            currMovie.windows.status = 0;
            continue;
        end
    end
    
    % STEP 5: Sample the protrusion vector in each window

    if ~isfield(currMovie,'protrusion') || ~isfield(currMovie.protrusion,'samples') || ...
            ~isfield(currMovie.protrusion.samples,'status') || ...
            currMovie.protrusion.samples.status ~= 1 || forceRun(5)
        try
            disp(['Sampling protrusion in movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieProtrusionSamples(currMovie,['protSamples_' ...
                winMethod '_' windowString  '.mat'],10,100, batchMode);
            
            if isfield(currMovie.protrusion.samples,'error')
               currMovie.protrusion.samples = rmfield(currMovie.protrusion.samples,'error');
           end
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);           
            currMovie.protrusion.samples.error = errMess;
            currMovie.protrusion.samples.status = 0;
            continue;
        end
        
    end 

    % STEP 6: Label
    
    if ~isfield(currMovie,'labels') || ~isfield(currMovie.labels,'status') || ...
            currMovie.labels.status ~= 1 || forceRun(6)
        try
            disp(['Labeling windows in movie ' num2str(iMovie) ' of ' num2str(nMovies)]);            
            currMovie = getMovieLabels(currMovie, batchMode);
            
            if isfield(currMovie.labels,'error')
                currMovie.labels = rmfield(currMovie.labels,'error');
            end
            
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.labels.error = errMess;
            currMovie.labels.status = 0;
            continue;
        end
    end
    
    % STEP 7: compute matrix D for figure 1
    
    if ~isfield(currMovie, 'output') || ~isfield(currMovie.output, 'fig1') || ...
            ~isfield(currMovie.output.fig1, 'status') || ...
            currMovie.output.fig1.status ~= 1 || forceRun(7)
        try
            disp(['Build Figure 1 data for movie ' num2str(iMovie) ' of ' num2str(nMovies)]);            
            currMovie = computeFigure1(currMovie, batchMode);
            
            if isfield(currMovie.output.fig1, 'error')
                currMovie.output.fig1 = rmfield(currMovie.output.fig1, 'error');
            end
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.output.fig1.error = errMess;
            currMovie.output.fig1.status = 0;
            continue;
        end
    end
    
    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMov) ': ' errMess.message], mfileName());
    end
    
    movieData{iMovie} = currMovie;
    
    disp([movieName ': DONE']);
end