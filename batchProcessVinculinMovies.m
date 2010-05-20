function batchProcessVinculinMovies(rootDirectory, params)

% params is a structure containing the following fields:
% .pixelSize
% .timeInterval
% .runSteps       a boolean array of size equal to the number of processes
% .batchMode      trigger for graphical display

procNames = {...
    'contours',...
    'protrusion',...
    'windows',...
    'protrusionSamples',...
    'labels',...
    'detection',...
    'tracking'};

procLocs = {...
    'contours',...
    'protrusion',...
    'windows',...
    'protrusion.samples',...
    'labels',...
    'detection',...
    'tracking'};

procFuns = {...
    'getMovieContours',...
    'getMovieProtrusion',...
    'getMovieWindows',...
    'getMovieProtrusionSamples',...
    'getMovieLabels',...
    'getMovieDetection',...
    'getMovieTracking'};

assert(numel(procNames) == numel(procLocs));
assert(numel(procNames) == numel(procFuns));

nProcesses = length(procFuns);

checkProcess = @(procName, movieData) ...
    eval(['checkMovie' upper(procName(1)) procName(2:end) '(movieData)']);

%
% CHECK INPUT ARGUMENTS
%

if nargin < 1 || isempty(rootDirectory)
    dataDirectory = uigetdir('', 'Select a data directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

if nargin < 2 || isempty(params)
    error('''params'' argument is missing.');
end

if ~isfield(params,'runSteps') || isempty(params.runSteps)
    params.runSteps = zeros(nSteps, 1);
end

if length(params.runSteps) ~= nProcesses
    error('size of parameter ''runSteps'' differs from number of processes.');
end

if ~isfield(params,'batchMode') || isempty(params.batchMode)
    params.batchMode = 1;
end

% Get every path from rootDirectory containing ch488 & ch560 subfolders.
paths = getDirectories(rootDirectory, 2, {'ch488', 'ch560'});

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

%
% BATCH ALL MOVIES
%

disp('Process all directories...');

nMovies = numel(paths);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    path = paths{iMovie};
   
    currMovie = movieData{iMovie};
    
    % INIT: SETUP MOVIE DATA

    try
        fieldNames = {...
            'bgDirectory',...
            'roiDirectory',...
            'tifDirectory',...
            'stkDirectory',...
            'analysisDirectory'};
        
        subDirNames = {'bg', 'roi', 'tif', 'stk', 'analysis'};
        
        channels = cell(numel(fieldNames), 1, 2);
        channels(:, 1, 1) = cellfun(@(x) [path filesep 'ch488' filesep x],...
            subDirNames, 'UniformOutput', false);
        channels(:, 1, 2) = cellfun(@(x) [path filesep 'ch560' filesep x],...
            subDirNames, 'UniformOutput', false);
        currMovie.channels = cell2struct(channels, fieldNames, 1);
        
        % We put every subsequent analysis in the ch488 analysis directory.
        currMovie.analysisDirectory = currMovie.channels(1).analysisDirectory;
        
        % Add these 2 fields to be compliant with Hunter's check routines:
        currMovie.imageDirectory = currMovie.channels(1).roiDirectory;
        currMovie.channelDirectory = {''};
        
        % Get the number of images
        
        n1 = numel(dir([currMovie.channels(1).roiDirectory filesep '*.tif']));
        n2 = numel(dir([currMovie.channels(2).roiDirectory filesep '*.tif']));
        
        % In case one of the channel hasn't been set, we still might want
        % to compute some processes.
        currMovie.nImages = max(n1,n2);
        
        assert(currMovie.nImages ~= 0);

        % Load physical parameter from
        filename = [currMovie.channels(2).analysisDirectory filesep 'fsmPhysiParam.mat'];
        
        if exist(filename, 'file')
            load(filename);

            if fsmPhysiParam.pixelSize ~= params.pixelSize
                disp([movieName ': pixel size in qFSM differs from value in batch params (SKIPPING).']);
                continue;
            end
        
            if fsmPhysiParam.frameInterval ~= params.timeInterval
                displ([movieName ': time interval in qFSM differs from value in batch params (SKIPPING).']);
                continue;
            end
            
            clear fsmPhysiParam;
        end
        
        currMovie.pixelSize_nm = params.pixelSize;
        currMovie.timeInterval_s = params.timeInterval;
        
        % Get the mask directory
        currMovie.masks.channelDirectory = {''};
        currMovie.masks.directory = [currMovie.channels(2).analysisDirectory...
            filesep 'edge' filesep 'cell_mask'];
        if exist(currMovie.masks.directory, 'dir')
            currMovie.masks.n = numel(dir([currMovie.masks.directory filesep '*.tif']));
            currMovie.masks.status = 1;
        else
            currMovie.masks.status = 0;
        end
        
        % Update from already saved movieData
        filename = [currMovie.analysisDirectory filesep 'movieData.mat'];
        if exist(filename, 'file')
            currMovie = load(filename);
            currMovie = currMovie.movieData;
        end
         
    catch errMess
        disp([movieName ': ' errMess.stack(1).name ':' num2str(errMess.stack(1).line) ' : ' errMess.message]);
        disp(['Error in movie ' num2str(iMovie) ': ' errMess.message '(SKIPPING)']);
        continue;
    end
    
    %
    % RUN PROCESSES
    %

    for iProc = 1:nProcesses
        procName = procNames{iProc};
        procFun = procFuns{iProc}; %#ok<NASGU>
        procLoc = procLocs{iProc};
        
        if  params.runSteps(iProc) == 1 || (~checkProcess(procName, currMovie) && params.runSteps(iProc) == 0)
            
            disp([movieName ': running process ' procName]);
            
            procParams = eval(['params.' procLoc]);
            procParamKeys = fieldnames(procParams);
            procParamsStr = arrayfun(@(iKey) ['procParams.' procParamKeys{iKey} ','], ...
                1:numel(procParamKeys), 'UniformOutput', false);
            % concatenate all params
            procParamsStr = strcat(procParamsStr{:}, 'params.batchMode');
            
            try
                currMovie = eval(['feval(procFun,currMovie,' procParamsStr ');']);
                
                if isfield(currMovie.(procName),'error')
                    eval(['currMovie.' procLoc ' = rmfield(currMovie.' procLoc ',''error'');']);
                end
 
            catch errMess
                disp(['Error in ' movieName ': ' errMess.stack(1).name ':' ...
                    num2str(errMess.stack(1).line) ' : ' errMess.message]);
                
                eval(['currMovie.' procLoc '.error = errMess;']);
                eval(['currMovie.' procLoc '.status = 0;']);
                continue;
            end
        else
            disp([movieName ': process ' procName ' not run']);
        end
    end   
    
    % Save results
    try
        %Save the updated movie data
        updateMovieData(currMovie)
    catch errMess
        errordlg(['Problem saving movie data in movie ' num2str(iMov) ': ' errMess.message], mfileName());
    end
    
    movieData{iMovie} = currMovie;

    disp([movieName ': DONE']);
end

%
% Create output directory for figures
%

% outputDirectory = fullfile(rootDirectory,'figures');
% if ~exist(outputDirectory, 'dir')
%     mkdir(rootDirectory, 'figures');
% end
% 
% % prefix the rootDirectory
% selectedPaths = paths(9:-1:8);
% 
% % suffix ch488/analysis
% selectedPaths = cellfun(@(subDir) fullfile(subDir,'ch488','analysis'),...
%     selectedPaths, 'UniformOutput', false);

%
% Make Figure 1
%
%disp('Make figure 1...');
%makeFAFigure1(selectedPaths, outputDirectory);

%
% Make Figure 2
%
%disp('Make figure 2...');
%makeFAFigure2(selectedPaths, outputDirectory);

%
% Make Figure 3
%
%disp('Make figure 3...');
%makeFAFigure3(selectedPaths, outputDirectory);


