function displayLayers(hFig, iFrame)

settings = get(hFig, 'UserData');

nLayers = numel(settings.layers);

% Get the axes of hFig.
hAxes = get(hFig, 'CurrentAxes');

% Get the channel plugins list
[dummy, layerPlugins] = getPlugins(); %#ok<ASGLU>

% Clear all layers
for iLayer = 1:nLayers
    tag = settings.layers{iLayer}.tag;
    
    hLayers = findall(get(hAxes, 'Children'), 'Tag', tag);
    
    delete(hLayers);
end

% Draw new layers
for iLayer = 1:nLayers
    layerTypeID = settings.layers{iLayer}.type;
    layerColor = settings.layers{iLayer}.color;
    tag = settings.layers{iLayer}.tag;
                
    % Display layer
    layerData = settings.layers{iLayer}.data{iFrame};
    
    layerPlugins(layerTypeID).displayFunc(hAxes, tag, layerData, layerColor);
end
