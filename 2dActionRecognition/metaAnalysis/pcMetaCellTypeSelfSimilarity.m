%% At the well level: EXPLAIN HERE
% D - similarity matrix
% ind - index of experiment (cell type) to be assessed
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function selfSimilarity = pcMetaCellTypeSelfSimilarity(D,cellTypeAll,ind)
    cellType = cellTypeAll{ind};    
    %     inds = strfind([cellTypeAll{:}],cellType); % WRONG!
    inds = find(strcmp(cellTypeAll,cellType));
    nInds = length(inds);
    
    selfSimilarity = [];
    
    for i = 1 : nInds
        if inds(i) == ind
            continue; % same experiment!
        end
        selfSimilarity = [selfSimilarity D(ind,inds(i))];
    end
end