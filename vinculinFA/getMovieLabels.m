function movieData = getMovieLabels(movieData, method, varargin)

validMethods = {'band', 'sector', 'window'};
isValidMethod = @(m) any(strcmp(m,validMethods));

ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('movieData', @checkMovieWindows);
ip.addRequired('method', isValidMethod);
ip.addParamValue('batchMode', true, @islogical);

ip.parse(movieData, method, varargin{:});
batchMode = ip.Results.batchMode;

%Indicate that labeling was started
movieData.labels.status = 0;

imSize = movieData.imSize;

%Load the windows
load([movieData.windows.directory filesep movieData.windows.fileName]);

%Check that they have been loaded
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

movieData.labels.dateTime = datestr(now);
movieData.labels.status = 1;

updateMovieData(movieData);

end