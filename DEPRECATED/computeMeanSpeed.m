function computeMeanSpeed
% COMPUTEMEANSPEED computes the mean speed within a region of interest. The
% ROI is defined by an image mask which has to be located in a sub
% directory called 'speedMask'. The chosen speedmap file to look at is
% defined by the number in the mask file name (e.g mask12.tif).

maskSpeedDir = uigetdir('', 'Select a mask directory:');

if ~ischar(maskSpeedDir)
    return;
end

maskFileNames = dir([maskSpeedDir filesep '*.tif']);

if ~numel(maskFileNames)
    error('Unable to locate mask files.');
end

speedMapDir = uigetdir('', 'Select the speed map directory:');

if ~ischar(speedMapDir)
    return;
end

speedMapFileNames = dir([speedMapDir filesep 'speedMap*.mat']);

if ~numel(speedMapFileNames)
    error('Unable to locate speed map files.');
end

[path, body, offset] = getFilenameBody(speedMapFileNames(1).name);
offset = -str2double(offset) + 1;

for iFile = 1:numel(maskFileNames)
    
    fileName = [maskSpeedDir filesep maskFileNames(iFile).name];
    BW = imread(fileName);
    
    [path, body, no] = getFilenameBody(fileName);    
    no = str2double(no);
    speedMapFileName = [speedMapDir filesep...
        speedMapFileNames(no + offset).name];
    load(speedMapFileName);
    
    if size(BW) ~= size(speedMap)
      error('Speed mask and speed map size differ.');
    end
    
    ind = find(BW);
    m = mean(speedMap(ind));
    disp([maskFileNames(iFile).name ': mean speed in mask = ' num2str(m)]);
    
end

end