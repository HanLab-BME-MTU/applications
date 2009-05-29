function displayLayers(hFig, iFrame)

settings = get(hFig, 'UserData');

numLayers = settings.numChannels;

% Get the axes of hFig.
hAxes = get(hFig, 'CurrentAxes');
hContent = get(hAxes, 'Children');

 % Should be read from a configuration file.
layerDisplayers = {@displaySpeckles};
layerTags = {'Cands'};

for iLayer = 1:numLayers
    layerType = settings.layers{iLayer}.type;
    layerColor = settings.layers{iLayer}.color;
    
    % Clear previous layer
    hLayers = findobj(hContent, '-regexp', 'Tag', layerTags{layerType});
    for h=1:length(hLayers)
        delete(hLayers(h));        
    end
                
    % Display layer
    fileName = [settings.layers{iLayer}.path ...
        filesep settings.layers{iLayer}.fileNames{iFrame}];
    
    layerDisplayers{layerType}(hAxes, fileName, layerColor);
end

end