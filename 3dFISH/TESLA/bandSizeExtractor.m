function bandSizeExtractor(bandStat, laneNum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xPos = zeros(numel(bandStat), 2);
size = zeros(numel(bandStat), 1);
for i = 1:numel(bandStat)
    xPos(i,1) = bandStat(i).bandPos(1);
    size(i,1) = bandStat(i).bandSize;
end

xPos(:,2) = kmeans(xPos(:, 1),laneNum);
indexNum = xPos(:,2);
bandHist = hist(indexNum, unique(indexNum));
bandSize = NaN(max(bandHist), laneNum);

count = 0;
for lane = 1:laneNum
    bandNum = bandHist(lane);
    for band = 1:bandNum
        count = count + 1;
        bandSize(band, lane) = size(count);
    end
end
save('bandSize.mat', 'bandSize')

end

