function [settings, status] = getFsmDataViewerSettings(handles)

% get the background index
h = findobj(handles, 'Tag', 'listboxBackground');
settings.backgroundIndex = get(h, 'Value');

% get the selected file name for the red channel
h = findobj(handles, 'Tag', 'editRedChannel');
settings.channels{1}.selectedFileName = get(h, 'String');

% get the selected file name for the green channel
h = findobj(handles, 'Tag', 'editGreenChannel');
settings.channels{2}.selectedFileName = get(h, 'String');

% get the selected file name for the blue channel
h = findobj(handles, 'Tag', 'editBlueChannel');
settings.channels{3}.selectedFileName = get(h, 'String');

% get the selected file name for the mask
h = findobj(handles, 'Tag', 'editMask');
settings.mask.selectedFileName = get(h, 'String');

% get the color and the selected file name for layers
for iLayer = 1:5
    h = findobj(handles, 'Tag', ['editLayer', num2str(iLayer)]);
    settings.layers{iLayer}.selectedFileName = get(h, 'String');
    h = findobj(handles, 'Tag', ['pushButtonColorLayer', num2str(iLayer)]);
    settings.layers{iLayer}.color = get(h, 'BackgroundColor');
end

for iChannel = 1:3
    [path, fileList, status] = getFileList(settings.channels{iChannel}.selectedFileName);    
    if ~status
        return;
    end
    settings.channels{iChannel}.path = path;
    settings.channels{iChannel}.fileList = fileList;
end

[path, fileList, status] = getFileList(settings.mask.selectedFileName);
if ~status
    return;
end
settings.mask.path = path;
settings.mask.fileList = fileList;

for iLayer = 1:5
    [path, fileList, status] = getFileList(settings.layers{iLayer}.selectedFileName);    
    if ~status
        return;
    end
    settings.layers{iLayer}.path = path;
    settings.layers{iLayer}.fileList = fileList;
end

% Check the number of files in every channel to be equal or null
firstNonZeroChannel = 0;
index = 1;
while index <=3
    if numel(settings.channels{index}.fileList)
        firstNonZeroChannel = index;
        break;
    end
    index = index + 1;
end
settings.firstNonZeroChannel = firstNonZeroChannel;

if index <= 3 % due to break
    for i=index+1:3
        numFiles = numel(settings.channels{i}.fileList);
        if numFiles && numFiles ~= numel(settings.channels{firstNonZeroChannel}.fileList)
            status = 0;
            errordlg('Number of files in background channels differ.');
            return;
        end
    end
end

numMaskFiles = numel(settings.mask.fileList);

% Check mask is provided in case no channel are set.
if ~firstNonZeroChannel && ~numMaskFiles
    status = 0;
    errordlg('Mask must be provided in case no channel has been set.');
    return;
end

% Check the number of mask files to matche the number of files in channels.
if firstNonZeroChannel && numMaskFiles && numMaskFiles ~= ...
        numel(settings.channels{firstNonZeroChannel}.fileList)
    status = 0;
    errordlg('Number of mask files does not match the number of files in channels.');
    return;
end

settings.numImages = numMaskFiles;

if firstNonZeroChannel
    settings.numImages = numel(settings.channels{firstNonZeroChannel}.fileList);
end

% Check the number of data files in layers to be lesser or equal to numFiles.
for iLayer = 1:5
    if numel(settings.layers{iLayer}.fileList) > settings.numImages
        status = 0;
        errordlg('Number of data files in layers cannot be greater than the number of files in channels.');
        return;        
    end
end
    
status = 1;

end


% 'getFileList' gets the list of files that have the same body name
% and are located in the same directory that 'fileName'.
% TODO: use this function everywhere in qfsm.

function [path, fileList, status] = getFileList(fileName)

if isempty(fileName)
    path = '';
    fileList = {};
    status = 1;
    return;
end

[path, body, no, ext] = getFilenameBody(fileName);
fileList = dir([path, filesep, body, '*', ext]);

% Check if there is any file.
if isempty(fileList)
    status = 0;
    errordlg(['No file name containing ' body ' can be found in ' path '.']);
    return;
end

% Rearrange files according to their number
fileNumbers = zeros(1, length(fileList));
for i = 1:length(fileList)
    [dummy1, dummy2, no] = getFilenameBody(fileList(i).name);
    fileNumbers(i) = str2double(no);
end

status = 1;

end