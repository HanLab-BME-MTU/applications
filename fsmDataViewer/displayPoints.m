function displayPoints(hAxes, tag, filename, layerColor)

data = load(filename);
if ~isstruct(data)
    error([filename, ' does not contain any data.']);
end
fieldNames = fieldnames(data);
for iField = 1:numel(fieldNames)
    P = data.(fieldNames{iField});
    
    if (numel(size(P)) == 2 && size(P, 2) == 2)
        line(P(:,2),P(:,1),'LineStyle', 'none',...
            'Marker', '.','Color',layerColor,'MarkerSize',8,...
            'Tag', tag, 'Parent', hAxes);
        return;
    end
end

error([filename, ' does not contain any data.']);