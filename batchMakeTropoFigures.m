function batchMakeTropoFigures(dataDirectory, analysisDirectory, forceRun, batchMode)

% STEP 1: Setup movie data
% STEP 2: Get movie contours            (depends on STEP 1)
% STEP 3: Compute protrusion            (depends on STEP 1)
% STEP 4: Get movie windows             (depends on STEP 2)
% STEP 5: Sample protrusion             (depends on STEP 3, 4)
% STEP 6: Get movie labels              (depends on STEP 4)
% STEP 7: Compute distance transform    (depends on STEP 1)
% STEP 8: Get movie density             (depends on STEP 1, 6)

nSteps = 8;

if nargin < 1 || isempty(dataDirectory)
    dataDirectory = uigetdir('', 'Select a data directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

if nargin < 2 || isempty(analysisDirectory)
    analysisDirectory = uigetdir('', 'Select an analysis directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

if nargin < 3 || isempty(forceRun)
    forceRun = zeros(nSteps, 1);
end

if length(forceRun) ~= nSteps
    error('Invalid number of steps in forceRun (3rd argument).');
end

if nargin < 4 || isempty(batchMode)
    batchMode = 1;
end

% Get every movies folder that contains 2 FSM subfolders either called:
% Actin GFP
% Actin TMx
% TMx TMy

isValidFSMProject = @(x) exist(fullfile(x, 'lastProjSettings.mat'),'file');
%dataPaths.ActinGFP = getDirectories(dataDirectory, 2, {'Actin', 'GFP'}, isValidFSMProject);
%dataPaths.ActinTM2 = getDirectories(dataDirectory, 2, {'Actin', 'TM2'}, isValidFSMProject);
%dataPaths.ActinTM4 = getDirectories(dataDirectory, 2, {'Actin', 'TM4'}, isValidFSMProject);
%dataPaths.ActinTM5 = getDirectories(dataDirectory, 2, {'Actin', 'TM5NM1'}, isValidFSMProject);
dataPaths.TM4TM2 = getDirectories(dataDirectory, 2, {'TM4', 'TM2'}, isValidFSMProject);
%dataPaths.TM5TM2 = getDirectories(dataDirectory, 2, {'TM5NM1', 'TM2'}, isValidFSMProject);
%dataPaths.TM5TM4 = getDirectories(dataDirectory, 2, {'TM5NM1', 'TM4'}, isValidFSMProject);

dataPathsFull = struct2cell(dataPaths);
dataPathsFull = vertcat(dataPathsFull{:});

disp('List of directories:');

for iMovie = 1:numel(dataPathsFull)
    disp([num2str(iMovie) ': ' dataPathsFull{iMovie}]);
end

% Create analysis directory for each data set

sstr = numel(dataDirectory);
analysisPaths = structfun(@(f) cellfun(@(x) fullfile(analysisDirectory,...
    x(sstr+1:end)), f, 'UniformOutput', false), dataPaths, ...
    'UniformOutput', false);

analysisPathsFull = struct2cell(analysisPaths);
analysisPathsFull = vertcat(analysisPathsFull{:});

for iMovie = 1:numel(analysisPathsFull)
    if ~exist(analysisPathsFull{iMovie}, 'dir');
        mkdir(analysisPathsFull{iMovie});
    end
end

%
% Process
%

disp('Process all directories...');

nMovies = numel(dataPathsFull);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(nMovies)];
    
    dataPath = dataPathsFull{iMovie};
    analysisPath = analysisPathsFull{iMovie};
   
    currMovie = movieData{iMovie};
    
    % STEP 1: Create movieData
    currMovie.analysisDirectory = analysisPath;
    content = dir(dataPath);
    content = {content.name};
    ind = cellfun(@(x) exist([dataPath filesep x filesep ...
        'lastProjSettings.mat'], 'file') ~= 0, content);
    if ~nnz(ind) || nnz(ind) ~= 2
        disp([movieName ': Unable to find the FSM directories. (SKIPPING).']);
        continue;
    end    
    currMovie.fsmDirectory = cellfun(@(x) [dataPath filesep x], content(ind), ...
        'UniformOutput', false);
    
    % Make sure TMs is the first directory in the list
    content = content(ind);
    if ~strcmpi(content{2}, 'actin') || ~strcmpi(content{2}, 'GFP')
        tmp = currMovie.fsmDirectory{1};
        currMovie.fsmDirectory{1} = currMovie.fsmDirectory{2};
        currMovie.fsmDirectory{2} = tmp;
    end
    
    % Add these 2 fields to be compliant with Hunter's check routines:
    % We choose arbitrary TM channel to be referred by 'imageDirectory'
    currMovie.imageDirectory = [currMovie.fsmDirectory{1} filesep 'crop'];
    currMovie.channelDirectory = {''};
    
    % Find the number of images
    currMovie.nImages = numel(dir([currMovie.imageDirectory filesep '*.tif']));
    
    % Load fsmPhysiParam.mat file. We search in both FSM directories. If
    % both directories contains that file, we choose the file in the first
    % directory.
    ind = arrayfun(@(x) exist([currMovie.fsmDirectory{x} filesep ...
        'fsmPhysiParam.mat'], 'file') ~= 0, 1:2);
    if any(ind)
        ind = find(ind, 1, 'first');
        load([currMovie.fsmDirectory{ind} filesep 'fsmPhysiParam.mat']);
    else
        disp([movieName ': Unable to locate fsmPhysiParam.mat (SKIPPING).']);
        continue;
    end        
    currMovie.pixelSize_nm = fsmPhysiParam.pixelSize;
    currMovie.timeInterval_s = fsmPhysiParam.frameInterval;
    clear fsmPhysiParam;

    % Set the mask directory (Actin)
    currMovie.masks.directory = [currMovie.fsmDirectory{2} filesep 'edge' filesep 'cell_mask'];
    
     % Add this field to be compliant with Hunter's check routines:
    currMovie.masks.channelDirectory = {''};
    currMovie.masks.n = numel(dir([currMovie.masks.directory filesep '*.tif']));
    currMovie.masks.status = 1;
    
    if exist([currMovie.analysisDirectory filesep 'movieData.mat'], 'file') && ~forceRun(1)
       currMovie = load([currMovie.analysisDirectory filesep 'movieData.mat']);
       currMovie = currMovie.movieData;
    end

    % STEP 2: Get contours
    dContour = 500 / currMovie.pixelSize_nm;
    
    if ~checkMovieContours(currMovie) || forceRun(2)
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

            handles.batch_processing = 1; %batchMode;
            handles.directory_name = [currMovie.masks.directory];
            handles.result_directory_name = currMovie.analysisDirectory;
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
                
                currMovie.protrusion.directory = handles.outdir;                
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
    winSize = 500 / currMovie.pixelSize_nm; % ~ 0.5um
    nBands = (5000 / (currMovie.pixelSize_nm * dContour)); % ~5 um depth
    iOuter = 2;
    iInner = 4;
    winMethod = 'p';            
    windowString = [num2str(dContour) 'by' num2str(winSize) 'pix_' ...
                num2str(iOuter) '_' num2str(iInner)];
            
    if ~checkMovieWindows(currMovie) || forceRun(4)
        try
            currMovie = setupMovieData(currMovie);

            disp(['Windowing movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieWindows(currMovie,winMethod,winSize,...
                nBands,iOuter,iInner,[],[],...
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
    if ~checkMovieProtrusionSamples(currMovie) || forceRun(5)
        try
            disp(['Sampling protrusion in movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieProtrusionSamples(currMovie,['protSamples_' ...
                winMethod '_' windowString  '.mat'], batchMode);
            
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

    % STEP 6: Activity Label (pause = 0, protrusion = 1, retraction = 2)
    if ~checkMovieLabels(currMovie) || forceRun(6)
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
    
    % STEP 7: Save Distance transform
    if ~checkMovieBWDist(currMovie) || forceRun(7)
        try
            disp(['Compute distance transform ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieBWDist(currMovie, batchMode);
            
            if isfield(currMovie.bwdist,'error')
                currMovie.bwdist = rmfield(currMovie.bwdist,'error');
            end
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.bwdist.error = errMess;
            currMovie.bwdist.status = 0;
            continue;
        end
    end
    
    % STEP 8: get Movie density map
    if ~checkMovieDensityMap(currMovie) || forceRun(8)
        try
            disp(['Compute Density map for movie ' num2str(iMovie) ' of ' num2str(nMovies)]);
            currMovie = getMovieDensityMap(currMovie, batchMode);
            
            if isfield(currMovie.density,'error')
                currMovie.density = rmfield(currMovie.density,'error');
            end
        catch errMess
            disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
            currMovie.density.error = errMess;
            currMovie.density.status = 0;
            continue;
        end
    end
    
    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMovie) ': ' errMess.message]);
    end
    
    movieData{iMovie} = currMovie;
    
    disp([movieName ': DONE']);
end

%
% Create output directory for figures
%

outputDirectory = [analysisDirectory filesep 'figures'];
if ~exist(outputDirectory, 'dir')
    mkdir(analysisDirectory, 'figures');
end

% Save the list of dataPathsFull in the figure in a text file format.
fid = fopen([outputDirectory filesep 'listFolders.txt'], 'w');
fprintf(fid, '%s\n%s\n%s\n', dataPathsFull{:});
fclose(fid);

% Figure 4
%disp('Make figure 4...');
%makeTropoFigure4(analysisPaths, outputDirectory);
% Figure 5
%disp('Make figure 5...');
%makeTropoFigure5(analysisPaths, outputDirectory);
% Figure 5bis
disp('Make figure 5bis...');
makeTropoFigure5bis(analysisPaths.TM4TM2, outputDirectory);
