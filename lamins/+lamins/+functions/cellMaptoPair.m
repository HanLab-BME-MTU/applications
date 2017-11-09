function pairs = cellMaptoPair(c)
% function pairs = cellMaptoPair(c)
% Convert a cell of length L into a N x 2 matrix map
% where the first column corresponds to the indices of the cell c
% and the second column refers to the values in c
    [idx,c] = cellfun(@(i,x) deal( ones(size(x(:)))*i, x(:) ),reshape(num2cell(1:numel(c)),size(c)),c,'UniformOutput',false);
    pairs = [vertcat(idx{:}) vertcat(c{:})];
end
    
