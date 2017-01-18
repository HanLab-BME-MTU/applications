function cellVideoClassification = analyzePlasticStates(matPosHist,matNegHist,vecBoundingBoxSize,vecFragmentVideoNumLabel,vecVideoClassLabel)
%Trains a classifier on blebbing, dendritic, background and any other
%available classes. Classifies known plastic cells with the classifier to
%track cell state transitions. Groups the fragments back into videos.

%some magic constants
PLASTIC_CLASS_NUM_LABEL = 3;
K_FEATS = 10;

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

%get the principal components of the features
[matPrincomp,~] = getPrincompMIP(matPosHist,matNegHist);
%Classifies the plastic fragments
classOut = classifyPlasticStates(matPrincomp,vecFragmentVideoClassLabel,K_FEATS,PLASTIC_CLASS_NUM_LABEL);

%groups the fragments into videos and stores the groups into a cell array
cellVideoClassification = cell(length(find(vecVideoClassLabel == PLASTIC_CLASS_NUM_LABEL)),1);
vecPlasticLabels = vecFragmentVideoNumLabel(vecFragmentVideoClassLabel == PLASTIC_CLASS_NUM_LABEL,:);
vecUniquePlasticLabels = unique(vecPlasticLabels);
for i = 1:length(find(vecVideoClassLabel == PLASTIC_CLASS_NUM_LABEL))
    cellVideoClassification{i} = classOut(vecPlasticLabels == vecUniquePlasticLabels(i),:);
end
end


function [classOut] = classifyPlasticStates(matPrincomp,vecFragmentVideoClassLabel,K_FEATS,PLASTIC_CLASS_NUM_LABEL)
%Defines and splits the feature vectors into plastic-class test sets and
%nonplastic-class training sets. Trains a decision tree classifier and
%classifies the heldout plastic video fragments. 

%selects only the K_FEATS number of principal components
matFeature = matPrincomp(:,1:K_FEATS);

%split into training and test classes
featuremat_train = matFeature(vecFragmentVideoClassLabel ~= PLASTIC_CLASS_NUM_LABEL,:);
featuremat_test = matFeature(vecFragmentVideoClassLabel == PLASTIC_CLASS_NUM_LABEL,:);
frameclass_train = vecFragmentVideoClassLabel(vecFragmentVideoClassLabel ~= PLASTIC_CLASS_NUM_LABEL,:);

%Trains a decision tree and uses it to classify the test set
myTree = fitctree(featuremat_train,frameclass_train);
classOut = predict(myTree,featuremat_test);
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