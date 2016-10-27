%% Input: cell type, link cell type --> source (comes from the meta data file)
function curCellSource = getCellSource(curCellType,allCellTypes)
ind = find(not(cellfun('isempty', strfind(lower(allCellTypes.ids),lower(curCellType)))));
curCellSource = allCellTypes.source{ind};
end