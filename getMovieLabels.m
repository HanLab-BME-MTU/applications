function movieData = getMovieLabels(movieData, method, batchMode)

%Indicate that labeling was started
movieData.labels.status = 0;

imSize = movieData.imSize';

%Verify that the windowing has been performed
if ~checkMovieWindows(movieData)
    error('Must window movie before labeling windows.');
end

%Load the windows
load([movieData.windows.directory filesep movieData.windows.fileName]);

%Check that they loaded
if ~exist('allWinPoly','var') || isempty(allWinPoly) %#ok<NODEF>
    error(['Problem loading windows from ' movieData.windows.directory ...
        filesep movieData.windows.fileName],mfilename);
end

movieData.labels.directory = [movieData.analysisDirectory filesep 'labels'];

if ~exist(movieData.labels.directory, 'dir')
    mkdir(movieData.labels.directory);
end

%Determine number of windows/bands
[nBands,nSectors,nFrames] = size(allWinPoly);

movieData.labels.nSectors = nSectors;
movieData.labels.nBands = nBands;
movieData.labels.nFrames = nFrames;
movieData.labels.method = method;

% BACKWARD COMPATIBILITY
if isfield(movieData.labels,'nWindows')
    movieData.labels = rmfield(movieData.labels,'nWindows');
end

%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

if ~batchMode
    h = waitbar(0,'Please wait, labeling window frames...');
end

for iFrame = 1:nFrames
    winPoly = allWinPoly(:,:,iFrame);
    
    labels = createLabelsFromWindows(winPoly, imSize, method);

    imwrite(uint16(labels), [movieData.labels.directory filesep 'labels_' ...
        num2str(iFrame,fString) '.tif'], 'Compression', 'lzw');
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame/nFrames,h)
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

if ~checkMovieProtrusionSamples(movieData)
    error('Must sample protrusion before labeling windows.');
end

% NOTE: Add this classification into protrusionSamples variable. This
% should add into getProtrusionSamples.m (Talk to Hunter). The following
% variables will be add to protrusionSample variable:
%
% stateNames
% states
% statePersistence

% BEGIN

% Classify edge velocity into 3 states: Pause, Protrusion, Retraction

%Load protrusion sample
load(fullfile(movieData.protrusion.directory, movieData.protrusion.samples.fileName));

%Check that the protrusion sample map loaded
if ~exist('protrusionSamples','var')
    error(['Problem loading protrusionSamples from ' movieData.windows.directory ...
        filesep movieData.windows.fileName],mfilename);
end

if ~isfield(protrusionSamples,'stateNames') || ~isfield(protrusionSamples,'states') || ...
        ~isfield(protrusionSamples,'statePersistence') %#ok<NODEF>
    
    %
    % Store state names
    %
    
    % Pause = 1, Protrusion = 2, Retraction = 3
    stateNames = {'Pause', 'Protrusion', 'Retraction'};
    protrusionSamples.stateNames = stateNames;
    
    %
    % classify speeds
    %
    
    states = ones(size(protrusionSamples.averageNormalComponent));
    states(protrusionSamples.averageNormalComponent >= 1) = 2;
    states(protrusionSamples.averageNormalComponent <= -1) = 3;
    
    protrusionSamples.states = states;
    
    %
    % compute state persistence
    %

    statePersistence = zeros(size(states));
    
    % Forward sweep: initialize 1st frame so that protrution / retraction
    % events which start at first frame are discarded.
    statePersistence(:,1) = 0;

    for iFrame = 2:nFrames-1
        for iState = 1:numel(stateNames)
            % iState at frames iFrame and iFrame-1
            ind = states(:,iFrame) == iState & ...
                states(:,iFrame-1) == iState & ...
                statePersistence(:,iFrame-1);
            statePersistence(ind,iFrame) = statePersistence(ind,iFrame-1) + 1;
    
            % iState at frame iFrame but not at iFrame-1
            ind = states(:,iFrame) == iState & states(:,iFrame-1) ~= iState;
            statePersistence(ind,iFrame) = 1;
        end
    end

    % Backward sweep: discard protrusion / retraction events which continue
    % at last frame.
    statePersistence(:,end) = 0;

    for iFrame = nFrames-2:-1:1
        % any state at frame iFrame and iFrame+1
        ind = states(:,iFrame) == states(:,iFrame+1) & ...
            ~statePersistence(:,iFrame+1);
        statePersistence(ind,iFrame) = 0;
    end
    
    protrusionSamples.statePersistence = statePersistence;
    
    % Update protrusionSamples file
    save(fullfile(movieData.protrusion.directory, ...
        movieData.protrusion.samples.fileName), 'protrusionSamples');
end

% END

movieData.labels.dateTime = datestr(now);
movieData.labels.status = 1;

updateMovieData(movieData);

end