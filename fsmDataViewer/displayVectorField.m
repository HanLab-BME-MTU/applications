function displayVectorField(hAxes, tag, filename, layerColor)

data = load(filename);
if ~isstruct(data)
    error([filename, ' does not contain any data.']);
end
fieldNames = fieldnames(data);
for iField = 1:numel(fieldNames)
    VF = data.(fieldNames{iField});
    
    if (numel(size(VF)) == 2 && size(VF, 2) == 4)
        set(hAxes, 'NextPlot', 'add');
        
        quiver(hAxes, VF(:, 2), VF(:, 1), VF(:, 4) - VF(:, 2), VF(:, 3) - VF(:, 1), 'Tag', tag, 'Color', layerColor);

        return;
    end
end

error([filename, ' does not contain any data.']);