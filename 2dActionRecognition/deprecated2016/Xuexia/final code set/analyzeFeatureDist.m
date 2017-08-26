function [matPosDist,matNegDist] = analyzeFeatureDist(matPosHist,matNegHist,vecBoundingBoxSize,vecFragmentVideoNumLabel)
%normalize feature histogram with regard to image size
for j = 1:length(vecBoundingBoxSize(:,1))
    matPosHist(j,:,:) = matPosHist(j,:,:)/vecBoundingBoxSize(j,1);
    matNegHist(j,:,:) = matNegHist(j,:,:)/vecBoundingBoxSize(j,1);
end

%split the classes
cellPosHistSplit = cell(length(unique(vecFragmentVideoNumLabel)),1);
cellNegHistSplit = cell(length(unique(vecFragmentVideoNumLabel)),1);
for l = 1:length(unique(vecFragmentVideoNumLabel))
    cellPosHistSplit{l} = matPosHist(squeeze(vecFragmentVideoNumLabel) == l ,:,:);
    cellNegHistSplit{l} = matNegHist(squeeze(vecFragmentVideoNumLabel) == l ,:,:);
end

%calculate centroids
pos_centroids = zeros(length(unique(vecFragmentVideoNumLabel)),length(matPosHist(1,:,1)),length(matPosHist(1,1,:)));
neg_centroids = zeros(length(unique(vecFragmentVideoNumLabel)),length(matNegHist(1,:,1)),length(matNegHist(1,1,:)));

for k = 1:length(unique(vecFragmentVideoNumLabel));
    pos_centroids(k,:,:) = mean(cellPosHistSplit{k});
    neg_centroids(k,:,:) = mean(cellNegHistSplit{k});
end

%calcualte distance between centroids
matPosDist = zeros(length(unique(vecFragmentVideoNumLabel)),length(unique(vecFragmentVideoNumLabel)),8);
matNegDist = zeros(length(unique(vecFragmentVideoNumLabel)),length(unique(vecFragmentVideoNumLabel)),8);
for m = 1:length(unique(vecFragmentVideoNumLabel))
    for n = 1:length(unique(vecFragmentVideoNumLabel))
        matPosDist(m,n,:) = sum(abs(pos_centroids(m,:,:) - pos_centroids(n,:,:)),3);
        matNegDist(m,n,:) = sum(abs(neg_centroids(m,:,:) - neg_centroids(n,:,:)),3);
    end
end
end