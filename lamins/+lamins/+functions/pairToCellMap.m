function cellMap = pairToCellMap(pairs)
%function cellMap = pairToCellMap(pairs)
% Converts a Nx2 pair matrix into a cell array where the indices correspond to
% unique values in the first column and 
% the corresponding cell refers to the corresponding values in the second column
    cellMap = accumarray(pairs(:,1),pairs(:,2),[],@(x) {x});
end

