%clc; clear all

import lamins.functions.*;
import lamins.classes.*;

dvpath1 = 'data/ag_072612_wt_Reconstructed 2.nd2';
dvpath2 = 'data/ag_080712wt_Reconstructed 3.nd2';

%%

MD = MovieData.load(dvpath1);
lamins = LaminsData(MD);
laminsCellRead = lamins.cellReader;
laminsCell = laminsCellRead.toCell;
zIdx = lamins.params.goodZ;
%lamins.showThumbs

laminsCellGood = squeeze(laminsCell(:,1,zIdx));
laminsCellGood([1 2],:) = laminsCellGood([2 1],:);

%% Make initial shadow segmentation

tempCell1 = cellfun(@(I) imadjust(I,stretchlim(I,0),[]), laminsCellGood, 'UniformOutput', 0);
tempCell1 = cellfun(@shadowSeg, tempCell1, 'UniformOutput', 0);

% tempCell = BW1 from the original implementation of nucleusSeg

% tempMat = cell2mat(tempCell1);
% figure; imshow(tempMat, 'InitialMagnification', 'fit')

% testIm = laminsCellGood{2,12};
% testIm = imadjust(testIm,stretchlim(testIm,0),[]);
% figure; imshow(imadjust(testIm), 'InitialMagnification', 'fit'); colormap jet 

% test1 = tempCell1{2,12};
% figure; imshow(test1, [], 'InitialMagnification', 'fit')

%% Close & fill

se = strel('disk', 10);
tempCell2 = cellfun(@(BW) imclose(BW,se), tempCell1, 'UniformOutput',0); %BW1c
tempCell2 = cellfun(@(BW) getForeground(~imfill(BW,'holes'),8), tempCell2, 'UniformOutput', 0); %BW1f

tempMat = cell2mat(tempCell2);
figure; imshow(tempMat, 'InitialMagnification', 'fit')


%% Create a boolean array that identifies good ellipses

setDim = size(tempCell1);
totNum = setDim(1)*setDim(2);
imDim = size(tempCell1{1,1});

sizes = cellfun(@(BW) sum(BW(:)), tempCell2);
[N, X] = histogram(sizes);
% figure; bar(X, N)

EVA = evalclusters(sizes(:),'kmeans','gap', 'KList', [1:3]);
opK = EVA.OptimalK;

[idx, cent] = kmeans(sizes(:), opK, 'Replicates', 10); 
[dummy, centIdx] = max(cent); % here we want the cluster with the largest mean
boolS1 = false(setDim); % first size clustered boolean
boolS1(idx==centIdx) = 1;

boolS1

%% Segment the easy ones

BWCell = cell(setDim);
for j=1:totNum    
    BWCell{j} = false(imDim);
end

% equivalent to BW1(~BWf) = 1
BWCell(boolS1) = cellfun(@(BW,BWf) BW|~BWf, tempCell1(boolS1), tempCell2(boolS1), 'UniformOutput', 0);
BWCell(boolS1) = cellfun(@(BW) imfill(getForeground(BW, 8), 'holes'), BWCell(boolS1), 'UniformOutput', 0);
BWCell(boolS1) = cellfun(@(BW) bwconvhull(BW), BWCell(boolS1), 'UniformOutput', 0);

% testIm = laminsCellGood{3,4};
% testIm = imadjust(testIm,stretchlim(testIm,0),[]);
% test1 = BWCell{3,4};
% figure; imshow(testIm, 'InitialMagnification', 'fit'); colormap jet
% hold on
% contour(test1, 'k', 'LineWidth', 2)

tempMat = cell2mat(BWCell);
figure; imshow(tempMat, 'InitialMagnification', 'fit')

%% Check if any easy ones failed

sizes = cellfun(@(BW) sum(BW(:)), BWCell(boolS1));

EVA = evalclusters(sizes(:),'kmeans','gap', 'KList', [1:3]);
opK = EVA.OptimalK;

[idx, cent] = kmeans(sizes(:), opK, 'Replicates', 10);
[dummy, centIdx] = max(cent);

listS1 = find(boolS1);
listS2 = listS1(idx~=centIdx);
boolS2 = false(setDim);
boolS2(listS2) = 1;

boolS2

%% Only for nonempty boolS2

BWCell(boolS2) = cellfun(@(BW,BWf) ~BW&BWf, tempCell1(boolS2), tempCell2(boolS2), 'UniformOutput',0);
BWCell(boolS2) = cellfun(@(BW) imopen(BW, se), BWCell(boolS2), 'UniformOutput', 0);
BWCell(boolS2) = cellfun(@(BW) bwconvhull(BW), BWCell(boolS2), 'UniformOutput', 0);

tempMat = cell2mat(BWCell);
figure; imshow(tempMat, 'InitialMagnification', 'fit')

%% Calculate max projection of the hard ones

for j=1:totNum % just to make sure it works    
    tempCell2{j} = false(imDim);
end

nBoolS1 = ~boolS1;
[xNS1, yNS1] = find(nBoolS1);

% If the entire column is empty this will do nothing
[listY, ia, ic] = unique(yNS1);

for j = 1:length(listY); % fill empty cells with the max projection from that column
    k = listY(j);
    boolRow = ic==j;
    maxProj = any(cat(3,BWCell{:,k}),3);
    tempCell2(xNS1(boolRow),k) = repmat({maxProj}, sum(boolRow), 1);
end

% tempMat = cell2mat(tempCell2);
% figure; imshow(tempMat, 'InitialMagnification', 'fit')


%% Segment the hard images

% test1 = tempCell1{1,1};
% figure; imshow(test1, [], 'InitialMagnification', 'fit')

BWCell(nBoolS1) = cellfun(@(BW,BWf) BW&BWf, tempCell1(nBoolS1), tempCell2(nBoolS1), 'UniformOutput',0);
BWCell(nBoolS1) = cellfun(@(BW) bwconvhull(BW), BWCell(nBoolS1), 'UniformOutput', 0);

% tempMat = cell2mat(BWCell);
% figure; imshow(tempMat, 'InitialMagnification', 'fit')

% testIm = laminsCellGood{1,1};
% testIm = imadjust(testIm,stretchlim(testIm,0),[]);
% test1 = BWCell{1,1};
% figure; imshow(testIm, 'InitialMagnification', 'fit'); colormap jet
% hold on
% contour(test1, 'k', 'LineWidth', 2)


%% Check to see if any hard ones failed

sizes = cellfun(@(BW) sum(BW(:)), BWCell(:));

EVA = evalclusters(sizes(:),'kmeans','gap', 'KList', [1:3]);
opK = EVA.OptimalK;

[idx, cent] = kmeans(sizes(:), opK, 'Replicates', 10);
[dummy, centIdx] = max(cent);

boolS3 = false(setDim);
boolS3(idx~=centIdx) =  1;
boolS3

%% Replace those with the max projection

BWCell(boolS3) = tempCell2(boolS3);
tempMat = cell2mat(BWCell);
figure; imshow(tempMat, 'InitialMagnification', 'fit')

%%

BWCell = cell(setDim);
for j=1:totNum    
    BWCell{j} = false(imDim);
end

BWCell(nBoolS1) = cellfun(@(BW,BWf) BW&BWf, tempCell1(nBoolS1), tempCell2(nBoolS1), 'UniformOutput',0);
BWCell(nBoolS1) = cellfun(@(BW) bwconvhull(BW), BWCell(nBoolS1), 'UniformOutput', 0);

tempMat = cell2mat(BWCell);
figure; imshow(tempMat, 'InitialMagnification', 'fit')