function cellArrayOut = LCH_filterEmptyCellArray(cellArrayIn)
cellArrayOut = {};
for i = 1 : length(cellArrayIn)
    if ~isempty(cellArrayIn{i})
        cellArrayOut{length(cellArrayOut) + 1} = cellArrayIn{i};
    end
end
end