function batchProcessMyMovies(params)

nProcesses = length(params.procNames);

checkFunc = @(procName, movieData) ...
    eval(['checkMovie' upper(procName(1)) procName(2:end) '(movieData)']);

procFunc = @(procName, paramFunc) ...
    eval(['getMovie' upper(procName(1)) procName(2:end) '(paramFunc{:})']);

%
% CHECK INPUT ARGUMENTS
%

if nargin < 1 || isempty(params.rootDirectory)
    dataDirectory = uigetdir('', 'Select a data directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

if nargin < 1 || isempty(params)
    error('''params'' argument is missing.');
end

if ~isfield(params,'runSteps') || isempty(params.runSteps)
    params.runSteps = zeros(nProcesses, 1);
end

if length(params.runSteps) ~= nProcesses
    error('size of parameter ''runSteps'' differs from number of processes.');
end

if ~isfield(params,'batchMode') || isempty(params.batchMode)
    params.batchMode = 1;
end

% Get every path from params.rootDirectory containing params.channelDirectory
% subfolders.
paths = getDirectories(params.rootDirectory, ...
    numel(params.channelDirectory), ...
    params.channelDirectory);

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

%
% BATCH ALL MOVIES
%

disp('Process all movies...');

nMovies = numel(paths);

movieData = cell(nMovies, 1);

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];
    
    path = paths{iMovie};
    
    % SETUP MOVIE DATA

    try
        currMovie = params.setupMovieDataFunc(path, params);
        
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
        procName = params.procNames{iProc};
        
        if  params.runSteps(iProc) == 1 || (~checkFunc(procName, currMovie) && params.runSteps(iProc) == 0)
            
            disp([movieName ': running process ' procName]);
            
            try
                procParams = struct2cell(eval(['params.' procName]));

                % save input parameters into the movieData
                eval(['currMovie.' procName '.params = params.' procName ';']);
                
                currMovie = procFunc(procName, vertcat({currMovie}, procParams, params.batchMode));

                if isfield(eval(['currMovie.' procName]),'error')                
                    eval(['currMovie.' procName ' = rmfield(currMovie.' procName ',''error'');']);
                end
 
            catch errMess
                disp(['Error in ' movieName ': ' errMess.stack(1).name ':' ...
                    num2str(errMess.stack(1).line) ' : ' errMess.message]);
                
                eval(['currMovie.' procName '.error = errMess;']);
                eval(['currMovie.' procName '.status = 0;']);
                continue;
            end
        else
            disp([movieName ': process ' procName ' not run']);
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
