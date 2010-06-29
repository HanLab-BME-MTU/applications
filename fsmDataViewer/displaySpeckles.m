function displaySpeckles(hAxes, tag, fileName, layerColor)

hImage = findobj(hAxes,'Type','image');
imSize = size(get(hImage,'CData'));

load(fileName);

if ~exist('cands', 'var')
    error('cands variable cannot be found in specified file.');
end

status = vertcat(cands(:).status); %#ok<NODEF>
cands = vertcat(cands(status == 1).Lmax);

% Get the list of speckles which had been already display
hSpeckles = findobj(hAxes,'-regexp','Tag','Speckles');

if ~isempty(hSpeckles)
    
    Z = false(imSize);
    ind = sub2ind(imSize, cands(:,1),cands(:,2));
    Z(ind) = true;
    
    interInd = cell(numel(hSpeckles),1);
    
    for iLayer = 1:numel(hSpeckles)
        xPos = get(hSpeckles(iLayer), 'XData');
        yPos = get(hSpeckles(iLayer), 'YData');
            
        ind2 = sub2ind(imSize, yPos, xPos);
        
        isSame = Z(ind2);
                
        % update current layer by removing intersecting speckles
        [Y,X] = ind2sub(imSize, ind2(~isSame));
        set(hSpeckles(iLayer), 'XData', X, 'YData', Y);
        
        interInd{iLayer} = ind2(isSame);
    end
    
    interInd = unique(horzcat(interInd{:}));
        
    % plot intersecting speckles using an yellow cross
    if any(interInd)
        [Y,X] = ind2sub(imSize, interInd);
        line(X,Y, 'LineStyle', 'none', 'Marker', 'x',...
            'Color',[1,1,0],'MarkerSize',8, 'Tag', tag, 'Parent', hAxes);
        
        % remove intersecting speckles from the list of cands
        Z(interInd) = false;
        [Y,X] = ind2sub(imSize,find(Z));
        cands = [Y,X];
    end
end

% plot the rest of the candidates
if ~isempty(cands)
    line(cands(:,2),cands(:,1),'LineStyle', 'none', 'Marker', '.',...
        'Color',layerColor,'MarkerSize',6, 'Tag', tag, 'Parent', hAxes);
end
