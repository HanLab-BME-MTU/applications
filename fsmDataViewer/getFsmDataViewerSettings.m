function [settings, status] = getFsmDataViewerSettings(handles)

status = 1;

%
% Get all data entered by the user
%

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
settings.mask.selectedFileNmae = get(h, 'String');

% get the color and the selected file name for layers
for iLayer = 1:5
    h = findobj(handles, 'Tag', ['editLayer', num2str(iLayer)]);
    settings.layers{iLayer}.selectedFileName = get(h, 'String');
    h = findobj(handles, 'Tag', ['pushButtonColorLayer', num2str(iLayer)]);
    settings.layers{iLayer}.color = get(h, 'BackgroundColor');
end

for iChannel = 1:3
    [fileList, status] = getFileList(settings.channels{iChannel}.selectedFileName);
    
    if ~status
        return;
    end
    
    settings.channels{iChannel}.fileList = fileList;
end

for iLayer = 1:5
    [fileList, status] = getFileList(settings.layers{iLayer}.selectedFileName);
    
    if ~status
        return;
    end

    settings.layers{iLayer}.fileList = fileList;
end

% Check that the number of files in every channel is the same
% (can be 0) TODO

%numFileChannel1 = numel(settings.channels{1}.fileList);
%numFileChannel2 = numel(settings.channels{2}.fileList);
%numFileChannel3 = numel(settings.channels{3}.fileList);
status = 0;
errordlg('Number of files in background channels differ.');
return;

end


% function get the list of files that have the same body name
% and are located in the same directory that 'fileName'.
% TODO: put it in common?

function [fileList status] = getFileList(fileName)

status = 1;

if ~isempty(fileName)
    fileList = {};
    return;
end

[path, body, no, ext] = getFilenameBody(fileName);
fileList = dir([path, filesep, body, '*.', ext]);

% Check if there is any file.
if isempty(fileList)
    status = 0;
    errordlg(['No file name containing ' bodyName ' can be found in ' path '.']);
    return;
end

% Rearrange files according to their number
fileNumbers = zeros(1, length(fileList));
for i = 1:length(fileList)
    [path, body, no] = getFilenameBody(fileList(i).name);
    fileNumbers(i) = str2double(no);
end

end