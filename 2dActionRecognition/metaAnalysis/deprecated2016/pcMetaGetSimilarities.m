
%% OBSOLITE!!
%% At the well level: EXPLAIN HERE
% D - similarity matrix
% ind - index of experiment (cell type) to be assessed
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function similarities = pcMetaGetSimilarities(D,inds)    
    nInds = length(inds);
    
    if nInds == 1
        similarities = nan;
        return;
    end
    
    similarities = nan(1,(nInds*(nInds-1))/2);
    
    ii = 1;
    for i = 1 : (nInds-1)
        for j = i+1 : nInds
            similarities(ii) = D(inds(i),inds(j));
            ii = ii + 1;
        end
    end
end