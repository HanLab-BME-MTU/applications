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
    [path, fileNames, status] = getFileNames(settings.channels{iChannel}.selectedFileName);    
    if ~status
        return;
    end
    settings.channels{iChannel}.path = path;
    settings.channels{iChannel}.fileNames = fileNames;
end

[path, fileNames, status] = getFileNames(settings.mask.selectedFileName);
if ~status
    return;
end
settings.mask.path = path;
settings.mask.fileNames = fileNames;

for iLayer = 1:5
    [path, fileNames, status] = getFileNames(settings.layers{iLayer}.selectedFileName);    
    if ~status
        return;
    end
    settings.layers{iLayer}.path = path;
    settings.layers{iLayer}.fileNames = fileNames;
end

% firstNonEmptyChannel corresponds to the index of the first non
% empty channel. May be 0 is not channel is provided.
settings.firstNonEmptyChannel = 0;
for iChannel = 1:3
    if numel(settings.channels{iChannel}.fileNames)
        settings.firstNonEmptyChannel = iChannel;
        break;
    end
end

if settings.firstNonEmptyChannel
    settings.numBackgroundFiles = ...
        numel(settings.channels{settings.firstNonEmptyChannel}.fileNames);
else
    settnigs.numBackgroundFiles = 0;
end

% Check the number of files in every channel to be equal or null
if settings.firstNonEmptyChannel
    for i=settings.firstNonEmptyChannel+1:3
        numFiles = numel(settings.channels{i}.fileNames);
        if numFiles && numFiles ~= settings.numBackgroundFiles
            status = 0;
            errordlg('Number of files in background channels differ.');
            return;
        end
    end
end

% Get the number of mask provided
settings.numMaskFiles = numel(settings.mask.fileNames);

% Check mask is provided in case no channel are set.
if ~settings.numBackgroundFiles && ~settings.numMaskFiles
    status = 0;
    errordlg('Mask must be provided in case no channel has been set.');
    return;
end

% Check every mask file is available for every background file.
if settings.numBackgroundFiles && settings.numMaskFiles
    for iFile = 1:settings.numBackgroundFiles
        [dummy, body, no] = ...
            getFilenameBody(settings.channels{settings.firstNonEmptyChannel}.fileNames{iFile});
        
        [dummy, found] = findNumberedFileInList(settings.mask.fileNames, str2double(no));

        if ~found
            status = 0;
            errordlg(['Unable to find the mask file for '...
                settings.channels{settings.firstNonEmptyChannel}.fileNames{iFile}, '.']);
            return;
        end
    end
end

% firstNonEmptyLayer corresponds to the index of the first non
% empty layer. May be 0 is not channel is provided.
settings.firstNonEmptyLayer = 0;
for iLayer = 1:5
    if numel(settings.layers{iLayer}.fileNames)
        settings.firstNonEmptyLayer = iLayer;
        break;
    end
end

if settings.firstNonEmptyLayer
    settings.numLayerFiles = ...
        numel(settings.layers{settings.firstNonEmptyLayer}.fileNames);
else
    settings.numLayerFiles = 0;
end

% Check the number of files in every layer to be equal or null
if settings.firstNonEmptyLayer
    for i=settings.firstNonEmptyLayer+1:3
        numFiles = numel(settings.layers{i}.fileNames);
        if numFiles && numFiles ~= settings.numLayerFiles
            status = 0;
            errordlg('Number of files in layers differ.');
            return;
        end
    end
end

% Check every mask file is available for every layer file.
if settings.numLayerFiles && settings.numMaskFiles
    for iFile = 1:settings.numLayerFiles
        [dummy, body, no] = ...
            getFilenameBody(settings.layers{settings.firstNonEmptyLayer}.fileNames{iFile});
        
        [dummy, found] = findNumberedFileInList(settings.mask.fileNames, str2double(no));

        if ~found
            status = 0;
            errordlg(['Unable to find the mask file for ' ...
                settings.layers{settings.firstNonEmptyLayer}.fileNames{iFile}, '.']);
            return;
        end
    end
end

% Check every layer file is available for every background file.
if settings.numBackgroundFiles && settings.numLayerFiles
    for iFile = 1:settings.numBackgroundFiles
        [dummy, body, no] = ...
            getFilenameBody(settings.channels{settings.firstNonEmptyChannel}.fileNames{iFile});
        
        no = str2double(no);
        
        [dummy, found] = findNumberedFileInList(settings.layers{firstNonEmptyLayer}.fileNames, str2double(no));
        
        if ~found
            status = 0;
            errordlg(['Unable to find the layer file for ' ...
                settings.channels{settings.firstNonEmptyChannel}.fileNames{iFile}, '.']);
            return;
        end
    end    
end

% Set the number of range.
settings.numFrames = settings.numBackgroundFiles;
if ~settings.numFrames
    settings.numFrames = settings.numLayerFiles;
    if ~settings.numFrames
        settings.numFrames = settings.numMaskFiles;
    end
end

status = 1;

end