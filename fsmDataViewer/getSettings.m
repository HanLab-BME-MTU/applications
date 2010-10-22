function [settings, status] = getSettings(hFig)

% Get names of channel type.
h = findobj(hFig, 'Tag', 'listboxChannelType');
channelTypeNames = get(h, 'String');

% Get channels
h = findobj(hFig, 'Tag', 'uitableChannels');
data = get(h, 'Data');
columnFormat = get(h, 'ColumnFormat');
colorNames = columnFormat{3};

nChannels = size(data, 1);

settings.channels = cell(nChannels,1);

for iChannel = 1:nChannels
    selected = data{iChannel, 1};
    
    if selected
        channelTypeName = data{iChannel, 2};
        channel.type = strcmp(channelTypeName, channelTypeNames) - 1;

        colorName = data{iChannel, 3};
        channel.color = strcmp(colorName, colorNames);        
        
        [path, fileNames] = getFileNames(data{iChannel, 4});
        
        channel.path = path;
        channel.fileNames = fileNames;
        
        settings.channels{iChannel} = channel;
    end
end

% update the number of channels
settings.channels = settings.channels(~cellfun(@isempty,settings.channels));
nChannels = numel(settings.channels);

% Check number of channels <= 3
if nChannels > 3
    status = 0;
    errordlg('Number of channels must be lesser or equal to 3.');
    return;
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

nLayers = size(data,1);

settings.layers = cell(nLayers,1);

for iLayer = 1:nLayers
    selected = data{iLayer, 1};
    
    if selected
        layerTypeName = data{iLayer, 2};
        layer.type = strmcmp(layerTypeName, layerTypeNames, 'exact') - 1;
        layer.tag = [layerTypeName '_' num2str(iLayer)];
        layer.color = data{iLayer, 3};
        
        [path, fileNames] = getFileNames(data{iLayer, 4});
        
        layer.path = path;
        layer.fileNames = fileNames;
        
        settings.layers{iLayer} = layer;
    end
end

% Update the number of layers
settings.layers = settings.layers(~cellfun(@isempty,settings.layers));
nLayers = numel(settings.layers);

% Check channel color compatibility
%
% - no multiple occurences of any color allowed

numColorsOcc = zeros(1, 3);
for iChannel = 1:settings.numChannels
    numColorsOcc(settings.channels{iChannel}.color) = ...
        numColorsOcc(settings.channels{iChannel}.color) + 1;
end

if any(numColorsOcc > 1)
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

% Get the number of frames (defined either by the number of image files or,
% if there is no channel provided, it is defined by the number of mask
% files.
if nChannels
    settings.nFrames = numel(settings.channels{1}.fileNames);
    
    nMasks = numel(settings.maskFileNames);
    
    if nMasks && nMasks ~= settings.nFrames
      status = 0;
      errordlg('Number of masks files differ from number of image files.');
      return
    end
    
elseif numel(settings.maskFileNames)
    settings.nFrames = numel(settings.maskFileNames);
else
  status = 0;
  errordlg('Mask must be provided in case no channel has been set.');
  return;
end

% Check the number of frames in every channel is equal
for iChannel = 1:nChannels
  if numel(settings.channels{iChannel}.fileNames) ~= settings.nFrames
    status = 0;
    errordlg('Number of files in channels differ.');
    return;
  end
end

% Check the number of frames in every layer is equal
for iLayer = 1:nLayers
  if numel(settings.layers{iLayer}.fileName) ~= settings.nFrames
    status = 0;
    errordlg('Number of files in layer differ.');
    return;
  end
end

status = 1;

end