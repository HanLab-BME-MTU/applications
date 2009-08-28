function displayLayers(hFig, iFrame)

settings = get(hFig, 'UserData');

numLayers = settings.numLayers;

% Get the axes of hFig.
hAxes = get(hFig, 'CurrentAxes');

% Get the channel plugins list
[channelPlugins layerPlugins] = getPlugins();

% Clear all layers
for iLayer = 1:numLayers
    tag = settings.layers{iLayer}.tag;
    
    hLayers = findall(get(hAxes, 'Children'), 'Tag', tag);
    
    delete(hLayers);
end

% Draw new layers
for iLayer = 1:numLayers
    layerTypeID = settings.layers{iLayer}.type;
    layerColor = settings.layers{iLayer}.color;
    tag = settings.layers{iLayer}.tag;
                
    % Display layer
    fileName = [settings.layers{iLayer}.path ...
        filesep settings.layers{iLayer}.fileNames{iFrame}];
    
    layerPlugins(layerTypeID).displayFunc(hAxes, tag, fileName, layerColor);
end

end