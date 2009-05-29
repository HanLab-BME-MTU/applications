function displayLayers(hFig, iFrame)

settings = get(hFig, 'UserData');

numLayers = settings.numChannels;

% Get the axes of hFig.
hAxes = get(hFig, 'CurrentAxes');
hContent = get(hAxes, 'Children');

% Get the channel plugins list
[channelPlugins layerPlugins] = getPlugins();

for iLayer = 1:numLayers
    layerTypeID = settings.layers{iLayer}.type;
    layerColor = settings.layers{iLayer}.color;
    tag = layerPlugins(layerTypeID).desc;
    
    % Clear previous layer
    hLayers = findobj(hContent, 'Tag', tag);
    for h=1:length(hLayers)
        delete(hLayers(h));        
    end
                
    % Display layer
    fileName = [settings.layers{iLayer}.path ...
        filesep settings.layers{iLayer}.fileNames{iFrame}];
    
    layerPlugins(layerTypeID).display(hAxes, tag, fileName, layerColor);
end

end