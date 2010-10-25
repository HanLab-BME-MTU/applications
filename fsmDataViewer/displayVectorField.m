function displayVectorField(hAxes, tag, VF, layerColor)

if (numel(size(VF)) == 2 && size(VF, 2) == 4)
    set(hAxes, 'NextPlot', 'add');
    
    quiver(hAxes, VF(:, 2), VF(:, 1), VF(:, 4), VF(:, 3), 0, 'Tag', tag, 'Color', layerColor);
        
    return;
end
