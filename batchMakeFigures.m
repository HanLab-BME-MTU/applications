function batchMakeFigures(varargin)

if nargin >= 1 && ~isempty(varargin{1})
    rootDirectory = varargin{1};
else
    % Ask for the root directory (Analysis2)
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

if nargin >= 2 && ~isempty(varargin{2})
    forceRun = varargin{2};
else
    forceRun = zeros(6, 1);
end

if nargin >= 3 && ~isempty(varargin{3})
    batchMode = varargin{3};
else
    batchMode = 1;
end

subFolders = {...
    ['TM2_Actin' filesep 'cell5'],...
    ['TM4_Actin' filesep '14_August_2009' filesep 'cell1'],...
    ['TM5NM1_Actin' filesep '26June2009' filesep 'Cell6'],...
    'TM2_TM4',...
    'TM2_TM5NM1',...
    'TM4_TM5NM1'};

paths = cellfun(@(x) [rootDirectory filesep x], subFolders,...
    'UniformOutput', false);

% not every subfolder are present now. Workaround:
ind = cellfun(@(x) exist(x, 'dir') ~= 0, paths);
paths = paths(ind);

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

disp('Process all directories...');

nMovies = numel(paths);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    path = paths{iMovie};
   
    currMovie = movieData{iMovie};
    
    % STEP 1: Create movieData
    currMovie.analysisDirectory = [path filesep 'windowAnalysis'];
    content = dir(path);
    content = {content.name};
    ind = cellfun(@(x) exist([path filesep x filesep ...
        'lastProjSettings.mat'], 'file') ~= 0, content);
    if ~nnz(ind) || nnz(ind) ~= 2
        disp([movieName ': Unable to find the FSM directories. (SKIPPING).']);
        continue;
    end    
    currMovie.fsmDirectory = cellfun(@(x) [path filesep x], content(ind), ...
        'UniformOutput', false);
    % Add these 2 fields to be compliant with Hunter's check routines:
    currMovie.imageDirectory = [currMovie.fsmDirectory{1} filesep 'crop'];
    currMovie.channelDirectory = {''};
    currMovie.nImages = numel(dir([currMovie.imageDirectory filesep '*.tif']));
    if exist([currMovie.fsmDirectory{1} filesep 'fsmPhysiParam.mat'], 'file')
        load([currMovie.fsmDirectory{1} filesep 'fsmPhysiParam.mat']);
    elseif exist([currMovie.fsmDirectory{2} filesep 'fsmPhysiParam.mat'], 'file')
        load([currMovie.fsmDirectory{2} filesep 'fsmPhysiParam.mat']);
    else
        disp([movieName ': Unable to locate fsmPhysiParam.mat (SKIPPING).']);
    end        
    currMovie.pixelSize_nm = fsmPhysiParam.pixelSize;
    currMovie.timeInterval_s = fsmPhysiParam.frameInterval;
    clear fsmPhysiParam;

    currMovie.masks.directory = [currMovie.fsmDirectory{1} filesep 'edge' filesep 'cell_mask'];
     % Add this field to be compliant with Hunter's check routines:
    currMovie.masks.channelDirectory = {''};
    currMovie.masks.n = numel(dir([currMovie.masks.directory filesep '*.tif']));
    currMovie.masks.status = 1;
    
    if exist([currMovie.analysisDirectory filesep 'movieData.mat'], 'file') && ~forceRun(1)
       currMovie = load([currMovie.analysisDirectory filesep 'movieData.mat']);
       currMovie = currMovie.movieData;
       
       % Remove previous outputs
       if isfield(currMovie, 'output')
           currMovie = rmfield(currMovie, 'output');
           updateMovieData(currMovie);
       end
    end

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

            handles.batch_processing = 1; %batchMode;
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
    nBandsLimit = 5; % ~ 5 um depth
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
        
    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMovie) ': ' errMess.message]);
    end
    
    movieData{iMovie} = currMovie;
    
    disp([movieName ': DONE']);
end

disp('Build figures...');

% Figure 3
makeFigure3(paths, batchMode);