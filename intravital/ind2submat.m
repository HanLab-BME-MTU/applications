function [submat] = ind2submat( matSize, indList )
    
    cellSubInd = cell(1, numel(matSize));    
    [cellSubInd{:}] = ind2sub(matSize, indList(:));
    submat = cell2mat( cellSubInd );
    
end