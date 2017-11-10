function [vector1AvgX,vector1AvgY,vector2AvgX,vector2AvgY] = compareQuiverVector(vector1,vector2,blockSize,imageSize)
%COMPAREQUIVERVECTOR compares average vector components in a given block size between two populations
% ideally, each vector should be 4x number of vectors array
yMap = nan(imageSize);
xMap = nan(imageSize);
a = ones(1,imageSize/blockSize)*blockSize;
[vector1AvgX,vector1AvgY,vector2AvgX,vector2AvgY]= deal(ones(length(a)));
ind = sub2ind([imageSize,imageSize],round(vector1(4,:)),round(vector1(3,:)));

yMap(ind) = vector1(2,:);
xMap(ind) = vector1(1,:);

yMapCell1 = mat2cell(yMap,a,a);
xMapCell1  =mat2cell(xMap,a,a);

% Do the same thing to vector 2
yMap = nan(imageSize);
xMap = nan(imageSize);
ind = sub2ind([imageSize,imageSize],round(vector2(4,:)),round(vector2(3,:)));

yMap(ind) = vector2(2,:);
xMap(ind) = vector2(1,:);

yMapCell2 = mat2cell(yMap,a,a);
xMapCell2  =mat2cell(xMap,a,a);

%Get average components from both vector maps

for k = 1:(imageSize/blockSize)
    for l = 1:(imageSize/blockSize)
        vector1AvgX(k,l) = nanmean(nanmean(xMapCell1{k,l}));
        vector1AvgY(k,l) = nanmean(nanmean(yMapCell1{k,l}));
        
        vector2AvgX(k,l) = nanmean(nanmean(xMapCell2{k,l}));
        vector2AvgY(k,l) = nanmean(nanmean(yMapCell2{k,l}));     
        
    end
end

vector1AvgX(isnan(vector1AvgX)) = 0;
vector1AvgY(isnan(vector1AvgY)) = 0;

vector2AvgX(isnan(vector2AvgX)) = 0;
vector2AvgY(isnan(vector2AvgY)) = 0;
end