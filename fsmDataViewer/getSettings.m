function [settings, status] = getSettings(hFig)

% Get names of channel type.
h = findobj(hFig, 'Tag', 'listboxChannelType');
channelTypeNames = get(h, 'String');

% Get channels
h = findobj(hFig, 'Tag', 'uitableChannels');
data = get(h, 'Data');
columnFormat = get(h, 'ColumnFormat');
colorNames = columnFormat{3};

settings.channels = {};

for iChannel = 1:size(data, 1)
    selected = data{iChannel, 1};
    
    if selected
        channelTypeName = data{iChannel, 2};
        channel.type = strmatch(channelTypeName, channelTypeNames) - 1;

        colorName = data{iChannel, 3};
        channel.color = strmatch(colorName, colorNames);        
        
        [path, fileNames] = getFileNames(data{iChannel, 4});
        
        channel.path = path;
        channel.fileNames = fileNames;
        
        settings.channels = vertcat(settings.channels, channel);
    end
end

% Get mask
settings.maskPath = '';
settings.maskFileNames = {};
h = findobj(hFig, 'Tag', 'checkboxMask');
value = get(h, 'Value');
if value
    h = findobj(hFig, 'Tag', 'editMask');
    
    [path, fileNames] = getFileNames(get(h, 'String'));
        
    settings.maskPath = path;
    settings.maskFileNames = fileNames;
end

% Get layers type names
h = findobj(hFig, 'Tag', 'listboxLayerType');
layerTypeNames = get(h, 'String');

% Get layers
h = findobj(hFig, 'Tag', 'uitableLayers');
data = get(h, 'Data');

settings.layers = {};

for iLayer = 1:size(data, 1)
    selected = data{iLayer, 1};
    
    if selected
        layerTypeName = data{iLayer, 2};
        layer.type = strmatch(layerTypeName, layerTypeNames) - 1;
        layer.tag = [layerTypeName '_' num2str(iLayer)];
        layer.color = data{iLayer, 3};
        
        [path, fileNames] = getFileNames(data{iLayer, 4});
        
        layer.path = path;
        layer.fileNames = fileNames;
        
        settings.layers = vertcat(settings.layers, layer);
    end
end

% Get the number of channels
settings.numChannels = numel(settings.channels);

% Check number of channels <= 3
if settings.numChannels > 3
    status = 0;
    errordlg('Number of channels must be lesser or equal to 3.');
    return;
end

% Check channel color compatibility
%
% - no multiple occurences of any color allowed

numColorsOcc = zeros(1, 3);
for iChannel = 1:settings.numChannels
    numColorsOcc(settings.channels{iChannel}.color) = ...
        numColorsOcc(settings.channels{iChannel}.color) + 1;
end

if max(numColorsOcc) > 1
    status = 0;
    errordlg('Multiple occurences of the same channel color is not possible.');
    return;
end

% Check channel type compatibility
if settings.numChannels > 1
    hasRawImage = 0;
    hasOtherData = 0;
    
    for iChannel = 1:settings.numChannels
        if settings.channels{iChannel}.type == 1
            hasRawImage = 1; % raw image
        else
            hasOtherData = 1; % other data
        end
    end
    
    if hasRawImage && hasOtherData
        status = 0;
        errordlg('Raw images can only be merged with raw image type.');
        return;
    end
end

% Get the number of layers
settings.numLayers = numel(settings.layers);

% Get the number of channel files
if settings.numChannels
    settings.numChannelFiles = numel(settings.channels{1}.fileNames);
else
    settings.numChannelFiles = 0;
end

% Check the number of files in every channel to be equal
for iChannel = 1:settings.numChannels
    if numel(settings.channels{iChannel}.fileNames) ~= settings.numChannelFiles
        status = 0;
        errordlg('Number of files in channels differ.');
        return;
    end
end

% Get the number of mask files
settings.numMaskFiles = numel(settings.maskFileNames);

% Get the number of layer files

% we would like to load layers with different number of files. e.g. speckles
% set, which are defined for every frame, with the vector fields, which are
% defined only on a subset of the whole sequence, due to temporal averaging.
% Hence we set 'numLayerFiles' to be the smallest number of files among the
% layers. 'iNumLayerFiles' corresponds to the argmin.

if settings.numLayers
    minNum = +Inf;
    for iLayer = 1:settings.numLayers
        curNum = numel(settings.layers{iLayer}.fileNames);
        if curNum < minNum
            minNum = curNum;
            iMinNum = iLayer;
        end
    end
    
    settings.numLayerFiles = minNum;
    settings.iNumLayerFiles = iMinNum;
else
    settings.numLayerFiles = 0;
    settings.iNumLayerFiles = 0;
end
        
% Check mask is provided in case no channel are set.
if ~settings.numChannelFiles && ~settings.numMaskFiles
    status = 0;
    errordlg('Mask must be provided in case no channel has been set.');
    return;
end

% Check every mask file is available for every channel file.
if settings.numChannels && settings.numMaskFiles
    for iFile = 1:settings.numChannelFiles
        [dummy, body, no] = ...
            getFilenameBody(settings.channels{1}.fileNames{iFile});
        
        [dummy, found] = findNumberedFileInList(settings.maskFileNames, str2double(no));

        if ~found
            status = 0;
            errordlg(['Unable to find the mask file for '...
                settings.channels{1}.fileNames{iFile}, '.']);
            return;
        end
    end
end

% Check every mask file is available for every layer file.
if settings.numLayers && settings.numMaskFiles
    for iFile = 1:settings.numLayerFiles
        [dummy, body, no] = ...
            getFilenameBody(settings.layers{settings.iNumLayerFiles}.fileNames{iFile});
        
        [dummy, found] = findNumberedFileInList(settings.maskFileNames, str2double(no));

        if ~found
            status = 0;
            errordlg(['Unable to find the mask file for ' ...
                settings.layers{settings.iNumLayerFiles}.fileNames{iFile}, '.']);
            return;
        end
    end
end

% Check every channel file is available for every layer file.
if settings.numLayers && settings.numChannels
    for iFile = 1:settings.numLayerFiles
        [dummy, body, no] = ...
            getFilenameBody(settings.layers{settings.iNumLayerFiles}.fileNames{iFile});
        
        no = str2double(no);
        
        [dummy, found] = findNumberedFileInList(settings.channels{1}.fileNames, no);
        
        if ~found
            status = 0;
            errordlg(['Unable to find the channel file for ' ...
                settings.layers{settings.iNumLayerFiles}.fileNames{iFile}, '.']);
            return;
        end
    end    
end

% set the number of frames
settings.numFrames = settings.numLayerFiles;
if ~settings.numFrames
    settings.numFrames = settings.numChannelFiles;
    if ~settings.numFrames
        settings.numFrames = settings.numMaskFiles;
    end
end

status = 1;

end