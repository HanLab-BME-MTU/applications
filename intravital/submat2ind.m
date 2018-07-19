function [ ind ] = submat2ind(imsize, submat)

    cellSubInd = cell(1, numel(imsize));
    for i = 1:numel(imsize)
        cellSubInd{i} = submat(:,i);        
    end
    
    ind = sub2ind(imsize, cellSubInd{:});
    
end