function movieData = getMovieActinFlowField(varargin)

% Input must be:
% varargin{1}: the movie data
movieData = varargin{1};
% varargin{2}: the channel index where the FSM Center analysis can be found
iChannel = varargin{2};
% varargin{3}: batch mode
batchMode = varargin{end};

%Indicate that particleDetection was started
movieData.actinFlowField.status = 0;

movieData.actinFlowField.directory = fullfile(movieData.analysisDirectory, ...
    'actinFlowField');

% Create output directory
if ~exist(movieData.actinFlowField.directory, 'dir')
    mkdir(movieData.actinFlowField.directory);
end

% Look for the file Md_#1-#2_d0=#3.mat
path = fullfile(movieData.imageDirectory, ...
    strtok(movieData.channelDirectory{iChannel},filesep), ...
    'analysis', 'post');

[status, result] = system(['find ' path ' -name ''Md*.mat''']);

assert(~status && numel(regexp(result,'\n')) == 1);

% Extract info from filename
result = result(1:end-1); % remove return carriage
tokens = regexp(result,filesep);
filename = result(tokens(end)+1:end);
matchstr = regexp(filename, '[0-9]+', 'match');
iFirst = str2double(matchstr{1});
iLast = str2double(matchstr{2});

% Load  Md_#1-#2_d0=#3.mat
load(result);

imSize = movieData.imSize;
nFrames = movieData.nImages;
%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

if ~batchMode
    h = waitbar(0, 'Please wait, retrieving actin flow field');
end
    
for iFrame = iFirst:iLast
    V = Md(:,:,iFrame-iFirst+1); %#ok<NODEF>
    dX = V(:,4)-V(:,2);
    dY = V(:,3)-V(:,1);
    theta = atan2(dY,dX);
    
    % Reshape
    thetaMap=reshape(theta,length(unique(V(:,2))),length(unique(V(:,1))))';
    thetaMap=imresize(thetaMap,imSize,'bicubic'); %#ok<NASGU>
    
    % Save to disk
    outputFile = fullfile(movieData.actinFlowField.directory,...
        ['thetaMap_' num2str(iFrame,fString) '.mat']);    
    save(outputFile,'thetaMap');
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame/(iLast-iFirst+1), h)
    end
end

% copy the leading frames (1...iFirst-1)
filename = fullfile(movieData.actinFlowField.directory,...
    ['thetaMap_' num2str(iFirst,fString) '.mat']);
load(filename);

for iFrame = 1:iFirst-1
    % Save to disk
    outputFile = fullfile(movieData.actinFlowField.directory,...
        ['thetaMap_' num2str(iFrame,fString) '.mat']);    
    save(outputFile,'thetaMap');
end

% copye the trailing frames (iLast+1...nFrames)
filename = fullfile(movieData.actinFlowField.directory,...
    ['thetaMap_' num2str(iLast,fString) '.mat']);
load(filename);

for iFrame = iLast+1:nFrames
    % Save to disk
    outputFile = fullfile(movieData.actinFlowField.directory,...
        ['thetaMap_' num2str(iFrame,fString) '.mat']);    
    save(outputFile,'thetaMap');
end

if ~batchMode && ishandle(h)
    close(h);
end

movieData.actinFlowField.dateTime = datestr(now);
movieData.actinFlowField.status = 1;

updateMovieData(movieData);
