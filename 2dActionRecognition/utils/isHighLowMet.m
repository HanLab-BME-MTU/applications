%% get metastatic efficiency
% input: cell type, 
% output: 'High' / 'Low' / nan
function curCellSource = isHighLowMet(curCellType,allCellTypes)
cellSource = getCellSource(curCellType,allCellTypes);
ind = find(not(cellfun('isempty', strfind(lower(allCellTypes.ids),lower(curCellType)))));
curCellSource = allCellTypes.source{ind};
end