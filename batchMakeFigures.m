function batchMakeFigures(dataDirectory, analysisDirectory, forceRun, batchMode)

nSteps = 7;

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
    error('Invalid number of steps in forceRun (2nd argument).');
end

if nargin < 4 || isempty(batchMode)
    batchMode = 1;
end

% Stick on this order: TM2, TM4, TM5, TM2_TM4, TM2_TM5, TM4_TM5
subFolders = {...
    ['TM2_Actin' filesep 'cell1_3'],...
    ['TM4_Actin' filesep '14_August_2009' filesep 'cell2'],...
    ['TM5NM1_Actin' filesep '26June2009' filesep 'Cell3'],...
    'TM2_TM4',...
    'TM2_TM5NM1',...
    'TM4_TM5NM1'};

dataPaths = cellfun(@(x) [dataDirectory filesep x], subFolders,...
    'UniformOutput', false);

% We process only files that exist now.
ind = cellfun(@(x) exist(x, 'dir') ~= 0, dataPaths);
dataPaths = dataPaths(ind);

disp('List of directories:');

for iMovie = 1:numel(dataPaths)
    disp([num2str(iMovie) ': ' dataPaths{iMovie}]);
end

analysisPaths = cellfun(@(x) [analysisDirectory filesep x], ...
    subFolders(ind), 'UniformOutput', false);

for iMovie = 1:numel(analysisPaths)
    if ~exist(analysisPaths{iMovie}, 'dir');
        mkdir(analysisPaths{iMovie});
    end
end

disp('Process all directories...');

nMovies = numel(dataPaths);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(nMovies)];
    
    dataPath = dataPaths{iMovie};
    analysisPath = analysisPaths{iMovie};
   
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
    if ~strcmpi(content{2}, 'actin')
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
    dContour = 500 / currMovie.pixelSize_nm; % ~ 500nm
    
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
    % Note: the width dWin should be set so that the autocorrelation over
    % windows of Edge Velocity Map is maximized (it might yield a trivial
    % solution dWin -> 0).
    dWin = 3000 / currMovie.pixelSize_nm; % ~ 3um
    iStart = 2;
    iEnd = 4;
    winMethod = 'e';            
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

    % STEP 6: Label
    nBandsLimit = 10; % ~5 um depth
    if ~isfield(currMovie,'labels') || ~isfield(currMovie.labels,'status') || ...
            currMovie.labels.status ~= 1 || forceRun(6)
        try
            disp(['Labeling windows in movie ' num2str(iMovie) ' of ' num2str(nMovies)]);            
            currMovie = getMovieLabels(currMovie, nBandsLimit, batchMode);
            
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
    if ~isfield(currMovie,'bwdist') || ~isfield(currMovie.bwdist,'status') || ...
            currMovie.bwdist.status ~= 1 || forceRun(7)
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
    
    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMovie) ': ' errMess.message]);
    end
    
    movieData{iMovie} = currMovie;
    
    disp([movieName ': DONE']);
end

outputDirectory = [analysisDirectory filesep 'figures'];
if ~exist(outputDirectory, 'dir')
    mkdir(analysisDirectory, 'figures');
end

% Save the list of dataPaths in the figure in a text file format.
fid = fopen([outputDirectory filesep 'listFolders.txt'], 'w');
fprintf(fid, '%s\n%s\n%s\n', dataPaths{:});
fclose(fid);

% Figure 4
disp('Make figure 4...');
makeFigure4(analysisPaths, outputDirectory);
% Figure 5
%disp('Make figure 5...');
%makeFigure5(analysisPaths, outputDirectory);