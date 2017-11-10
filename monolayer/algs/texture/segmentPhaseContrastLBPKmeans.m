%% Input:
% I - input image
% patchSize - size of patches
% lbpMapping
% 
% Output: binary image at original resolution

function [BIN,result] = segmentPhaseContrastLBPKmeans(I,patchSize,lbpMapping)

if size(I,3) == 3
    I = rgb2gray(I);
end

if nargin < 3
    lbpMapping = getmapping(8,'riu2');
end

Ilbp = lbp(I,1,8,lbpMapping,'');
[result] = segmentLbpImage(Ilbp,patchSize);
bin = false(1,size(result.pacthHist,1));
[cellsLabel,labelDiff] = findCellLabel(result);
bin(result.kmeans.labels == cellsLabel) = true;
BIN  = reshape(bin,size(result.IHist,1),size(result.IHist,2));
BIN = imresize(BIN,size(I));

if ((labelDiff < 0.2) && (std(double(I(BIN))) < std(double(I(~BIN)))))
   BIN = ~BIN; 
end
end

function [cellsLabel,labelDiff] = findCellLabel(result)
label1 = result.kmeans.normDistances(result.kmeans.labels == 1);
label2 = result.kmeans.normDistances(result.kmeans.labels == 2);

labelDiff = abs(mean(label1) - mean(label2));
if mean(label1) < mean(label2)
    cellsLabel = 1;
else
    cellsLabel = 2;
end
end