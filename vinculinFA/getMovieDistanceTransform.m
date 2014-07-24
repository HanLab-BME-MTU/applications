function movieData = getMovieDistanceTransform(movieData, varargin)

%Verify that the masks have been created
checkMovieData = @(movieData) checkMovieMasks(movieData, 1);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', checkMovieData);
ip.addParamValue('batchMode', true, @islogical);

ip.parse(movieData, varargin{:});
batchMode = ip.Results.batchMode;

%Indicate that computing distance transform was started
movieData.distanceTransform.status = 0;

movieData.distanceTransform.directory = [movieData.analysisDirectory ...
    filesep 'distanceTransform'];

if ~exist(movieData.distanceTransform.directory, 'dir')
    mkdir(movieData.distanceTransform.directory);
end

nFrames = movieData.nImages;
%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

% Read the list of Actin masks
maskPath = movieData.masks.directory;
maskFiles = dir([maskPath filesep '*.tif']);

%Go through each frame and save the windows to a file
if ~batchMode
    h = waitbar(0,'Please wait, computing distance transform...');
end

for iFrame = 1:nFrames    
    BW = imread([maskPath filesep maskFiles(iFrame).name]);
    distToEdge = double(bwdist(1 - BW)); %#ok<NASGU>

    save([movieData.distanceTransform.directory filesep 'distanceTransform_' ...
        num2str(iFrame,fString) '.mat'], 'distToEdge');
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame/nFrames,h);
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

movieData.distanceTransform.dateTime = datestr(now);
movieData.distanceTransform.status = 1;

updateMovieData(movieData);

end