
%% OBSOLITE!!
%% At the well level: EXPLAIN HERE
% D - similarity matrix
% ind - index of experiment (cell type) to be assessed
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function similarities = pcMetaGetInterSimilarities(D,inds1,inds2)    
    nInds1 = length(inds1);
    nInds2 = length(inds2);
    
    if nInds1 == 1
        similarities = nan;
        return;
    end
    
    similarities = nan(1,nInds1*nInds2);
    
    ii = 1;
    for i = 1 : nInds1
        for j = 1 : nInds2
            similarities(ii) = D(inds1(i),inds2(j));
            ii = ii + 1;
        end
    end
end