function displayLayers(hFig, iFrame)

settings = get(hFig, 'UserData');

numLayers = settings.numLayers;

% Get the axes of hFig.
hAxes = get(hFig, 'CurrentAxes');

% Get the channel plugins list
[channelPlugins layerPlugins] = getPlugins();

for iLayer = 1:numLayers
    layerTypeID = settings.layers{iLayer}.type;
    layerColor = settings.layers{iLayer}.color;
    tag = settings.layers{iLayer}.tag;

    % Clear previous layer
    hLayers = findall(get(hAxes, 'Children'), 'Tag', tag);
    for h=1:length(hLayers)
        delete(hLayers(h));        
    end
                
    % Display layer
    fileName = [settings.layers{iLayer}.path ...
        filesep settings.layers{iLayer}.fileNames{iFrame}];
    
    layerPlugins(layerTypeID).displayFunc(hAxes, tag, fileName, layerColor);
end

end