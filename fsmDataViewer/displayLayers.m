function displayLayers(hFig, iFrame)

settings = get(hFig, 'UserData');

% Get the number of non-empty layers.
numNonEmptyLayers = 0;
for iLayer = 1:3
    path = settings.layers{iLayer}.path;
    if ~isempty(path)
        numNonEmptyLayers = numNonEmptyLayers + 1;
    end
end

if numNonEmptyLayers
    % Clear the previous layers
    hAxes = get(hFig, 'CurrentAxes');
    hContent = get(hAxes, 'Children');
    % Clear previous speckles
    hLayers = findobj(hContent, '-regexp', 'Tag', '.Cands');
    for h=1:length(hLayers)
        delete(hLayers(h));        
    end
    
    % Display all the layers
    for iLayer = 1:5
        path = settings.layers{iLayer}.path;
        
        if ~isempty(path)
            % Get the color associated with the layer
            color = settings.layers{iLayer}.color;

            % Open the data.
            fileName = settings.layers{iLayer}.fileNames{iFrame};
            data = load([path, filesep fileName]);

            % Get the type of layer (speckle, vector flow, etc.)
            %layerType = getLayerType(data);
            layerType = 'speckles';
            
            switch layerType
                case 'speckles',
                    cands = data.cands;
                    clear data;

                    % Extract speckle classes
                    primary = find([cands.status] == 1 & [cands.speckleType] == 1);
                    secondary = find([cands.status] == 1 & [cands.speckleType] == 2);
                    tertiary = find([cands.status] == 1 & [cands.speckleType] == 3);
                    higher = find([cands.status] == 1 & [cands.speckleType] > 3);

                    % Extract speckle positions
                    pPos=reshape([cands(primary).Lmax], 2, length([cands(primary).Lmax])/2)';
                    sPos=reshape([cands(secondary).Lmax], 2, length([cands(secondary).Lmax])/2)';
                    tPos=reshape([cands(tertiary).Lmax], 2, length([cands(tertiary).Lmax])/2)';
                    hPos=reshape([cands(higher).Lmax], 2, length([cands(higher).Lmax])/2)';

                    % Primary speckles
                    h1 = line(pPos(:,2),pPos(:,1),'LineStyle', 'none',...
                        'Marker', '.','Color',color,'MarkerSize',6);                    
                    set(h1, 'Tag', 'pCands');
                    
                    % Secondary peckles
                    h2 = line(sPos(:,2),sPos(:,1),'LineStyle', 'none',...
                        'Marker', '+','Color',color,'MarkerSize',4); 
                    set(h2, 'Tag', 'sCands');
                    
                    % Tertiary speckles
                    h3 = line(tPos(:,2),tPos(:,1),'LineStyle', 'none',...
                        'Marker', '^','Color',color,'MarkerSize',4);
                    set(h3, 'Tag', 'tCands');
                    
                    % Higher-order speckles
                    h4 = line(hPos(:,2),hPos(:,1),'LineStyle', 'none',...
                        'Marker', '*','Color',color,'MarkerSize',4);
                    set(h4, 'Tag', 'hCands');
                    
                otherwise,
                    error('Cannot identify the type of layer from data file.');
            end            
        end
    end
end

end