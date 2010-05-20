function batchProcessTropoMovies(dataDirectory, analysisDirectory, params)

% params is a structure containing the following fields:
% .pixelSize
% .timeInterval
% .forceRun       a boolean array of size equal to the number of processes
% .batchMode      trigger for graphical display

procNames = {...
    'contours',...
    'protrusion',...
    'windows',...
    'protrusionSamples',...
    'labels',...
    'bwdist',...
    'density'};

procLocs = {...
    'contours',...
    'protrusion',...
    'windows',...
    'protrusion.samples',...
    'labels',...
    'bwdist',...
    'density'};

procFuns = {...
    'getMovieContours',...
    'getMovieProtrusion',...
    'getMovieWindows',...
    'getMovieProtrusionSamples',...
    'getMovieLabels',...
    'getMovieBWDist',...
    'getMovieDensityMap'};

assert(numel(procNames) == numel(procLocs));
assert(numel(procNames) == numel(procFuns));

nProcesses = length(procFuns);

checkProcess = @(procName, movieData) ...
    eval(['checkMovie' upper(procName(1)) procName(2:end) '(movieData)']);

%
% CHECK INPUT ARGUMENT
%

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

if nargin < 3 || isempty(params)
    error('''params'' argument is missing.');
end

if ~isfield(params,'forceRun') || isempty(params.forceRun)
    params.forceRun = zeros(nSteps, 1);
end

if length(params.forceRun) ~= nProcesses
    error('size of parameter ''forceRun'' differs from number of processes.');
end

if ~isfield(params,'batchMode') || isempty(params.batchMode)
    params.batchMode = 1;
end

% Get every movies folder that contains 2 FSM subfolders either called:
% Actin GFP
% Actin TMx
% TMx TMy

isValidFSMProject = @(x) exist(fullfile(x, 'lastProjSettings.mat'),'file');
%dataPaths.ActinGFP = getDirectories(dataDirectory, 2, {'Actin', 'GFP'}, isValidFSMProject);
dataPaths.ActinTM2 = getDirectories(dataDirectory, 2, {'Actin', 'TM2'}, isValidFSMProject);
% dataPaths.ActinTM4 = getDirectories(dataDirectory, 2, {'Actin', 'TM4'}, isValidFSMProject);
% dataPaths.ActinTM5 = getDirectories(dataDirectory, 2, {'Actin', 'TM5NM1'}, isValidFSMProject);
% dataPaths.TM4TM2 = getDirectories(dataDirectory, 2, {'TM4', 'TM2'}, isValidFSMProject);
% dataPaths.TM5TM2 = getDirectories(dataDirectory, 2, {'TM5NM1', 'TM2'}, isValidFSMProject);
% dataPaths.TM5TM4 = getDirectories(dataDirectory, 2, {'TM5NM1', 'TM4'}, isValidFSMProject);

% Concatenate all data paths
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

% Concatenate all data paths
analysisPathsFull = struct2cell(analysisPaths);
analysisPathsFull = vertcat(analysisPathsFull{:});

for iMovie = 1:numel(analysisPathsFull)
    if ~exist(analysisPathsFull{iMovie}, 'dir');
        mkdir(analysisPathsFull{iMovie});
    end
end

%
% BATCH ALL MOVIES
%

disp('Process all directories...');

nMovies = numel(dataPathsFull);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(nMovies)];
    
    dataPath = dataPathsFull{iMovie};
    analysisPath = analysisPathsFull{iMovie};
   
    currMovie = movieData{iMovie};
    
    %
    % INIT: SETUP MOVIE DATA
    %
    
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
        
        if fsmPhysiParam.pixelSize ~= params.pixelSize
            disp([movieName ': pixel size in qFSM differs from value in batch params (SKIPPING).']);
            continue;
        end
        
        if fsmPhysiParam.frameInterval ~= params.timeInterval
            displ([movieName ': time interval in qFSM differs from value in batch params (SKIPPING).']);
            continue;
        end
        
        clear fsmPhysiParam;        
    else
        disp([movieName ': Unable to locate fsmPhysiParam.mat (SKIPPING).']);
        continue;
    end
    
    currMovie.pixelSize_nm = params.pixelSize;
    currMovie.timeInterval_s = params.timeInterval;

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

    %
    % RUN PROCESSES
    %

    for iProc = 1:nProcesses
        procName = procNames{iProc};
        procFun = procFuns{iProc};
        procLoc = procLocs{iProc};
        
        if ~checkProcess(procName, currMovie) || params.forceRun(iProc)
            
            disp([movieName ': running process ' procName{iProc}]);
            
            procParams = eval(['params.' procLoc]);
            procParamKeys = fieldnames(procParams);
            procParamsStr = arrayfun(@(iKey) ['procParams.' procParamKeys{iKey} ','], ...
                1:numel(procParamKeys), 'UniformOutput', false);
            % concatenate all params
            procParamsStr = strcat(procParamsStr{:}, params.batchMode);
            
            try
                currMovie = procFun(currMovie, procParamsStr);
                
                if isfield(currMovie.(procName),'error')
                    eval(['currMovie.' procLoc ' = rmfield(currMovie.' procLoc ',''error'')']);
                end
 
            catch errMess
                disp(['Error in ' movieName ': ' errMess.stack(1).name ':' ...
                    num2str(errMess.stack(1).line) ' : ' errMess.message]);
                
                eval(['currMovie.' procLoc '.error = errMess']);
                eval(['currMovie.' procLoc '.status = 0']);
                continue;
            end
        else
            disp([movieName ': process ' procName{iProc} ' not run']);
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
%disp('Make figure 5bis...');
%makeTropoFigure5bis(analysisPaths.TM4TM2, outputDirectory);
