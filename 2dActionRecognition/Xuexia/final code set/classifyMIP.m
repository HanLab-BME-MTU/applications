function [cp1nn,cpDtree] = classifyMIP(matPosHist,matNegHist,vecBoundingBoxSize,vecFragmentVideoNumLabel,vecVideoClassLabel)
%normalize features respect to image size
for i = 1:length(vecBoundingBoxSize(:,1))
    matPosHist(i,:,:) = matPosHist(i,:,:)./vecBoundingBoxSize(i,1);
    matNegHist(i,:,:) = matNegHist(i,:,:)./vecBoundingBoxSize(i,1);
end

%set frame classes to video classes
vecFragmentVideoClassLabel = zeros(length(vecFragmentVideoNumLabel),1);
for j = 1:length(vecFragmentVideoClassLabel)
   vecFragmentVideoClassLabel(j) = vecVideoClassLabel(vecFragmentVideoNumLabel(j));
end

% performs principal component analysis for dimensionality reduction
matPrincomp = getPrincompMIP(matPosHist,matNegHist);

% cross validation
% all features,video labels, frame class, number of features to use, k-fold, number of repeats
[cp1nn,cpDtree] = getCrossvalAcc(matPrincomp,vecFragmentVideoNumLabel,vecFragmentVideoClassLabel,10,10,50);
end


function matPrincomp = getPrincompMIP(matPosHist,matNegHist)
%performs principal component analysis on the feature histograms for
%dimensonality reduction

%linearize the histogram matricies and concatenates the positive and
%negative MIP descriptor sets
cathistmatpos = squeeze(reshape(matPosHist,size(matPosHist,1),1,[]));
cathistmatneg = squeeze(reshape(matNegHist,size(matNegHist,1),1,[]));
cathistmat = [cathistmatpos,cathistmatneg];

[~,matPrincomp,~] = princomp(cathistmat);
end


function [cp1nn,cpDtree] = getCrossvalAcc(matPrincomp,vecFragmentVideoNumLabel,vecFragmentVideoClassLabel,K_FEATS,K_FOLD,N_REPEATS)
%performs K_FOLD cross validation using K_FEATS principal components 
%N_REPEATS times and returns a 1-nearest neighbor accuracy (cp1nn)
%and a decision tree accuracy (cpDtree).

%selects only the K_FEATS number of principal components
matFeatures = matPrincomp(:,1:K_FEATS);

%initialize the classperf objects
cp1nn = classperf(vecFragmentVideoClassLabel);
cpDtree = classperf(vecFragmentVideoClassLabel);

%peforms N_REPEATS replicates
for n = 1:N_REPEATS
    
    %randomly splits the cells into training and test groups
    vecKfoldCellSelection = crossvalind('Kfold',size(unique(vecFragmentVideoNumLabel),1),K_FOLD);
    for k = 1:K_FOLD
        %finds indicies of videos selected to be test set in current kth fold
        testidx = arrayfun(@(x) find(vecFragmentVideoNumLabel == x), find(vecKfoldCellSelection == k), 'UniformOutput', false);
        %converts from cell to matrix
        testidx = cell2mat(testidx);
        %transforms indicies to logical        
        test = zeros(size(vecFragmentVideoClassLabel));
        test(testidx) = 1;
        test = test > 0;
        %finds the training set
        train = ~test;
        
        %1nn classification
        class1nn = knnclassify(matFeatures(test,:),matFeatures(train,:),vecFragmentVideoClassLabel(train,:));
        classperf(cp1nn,class1nn,test);
        
        %decision tree classification
        mytree = fitctree(matFeatures(train,:),vecFragmentVideoClassLabel(train,:));
        classdtree = predict(mytree,matFeatures(test,:));
        classperf(cpDtree,classdtree,test);
    end
end
end