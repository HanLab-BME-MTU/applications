function [cellClassPrincomp] = analyzeFeaturePrincomp(matPosHist,matNegHist,vecBoundingBoxSize,vecFragmentVideoNumLabel,vecVideoClassLabel)
%performs principal component analysis for dimensionality reduction and
%splits the principal components into labeled classes

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

%performs principal component analysis for dimensionality reduction
matPrincomp = getPrincompMIP(matPosHist,matNegHist);

%split matPrincomp into cell types
cellClassPrincomp = cell(length(unique(vecVideoClassLabel)),1);
for k = 1:length(unique(vecVideoClassLabel))
    cellClassPrincomp{i} = matPrincomp(vecFragmentVideoClassLabel == k,:);
end
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
