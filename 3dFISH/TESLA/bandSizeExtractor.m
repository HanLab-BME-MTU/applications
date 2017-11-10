function bandSizeExtractor()
%BANDSIZEEXTRACTOR takes input TeSLA analysis folder and check for bandStat
%to generate records of band size for each lane
%   Detailed explanation goes here

[fileName, filePath] = uigetfile('*.mat', 'Select mat file with bandStat', ...
    'X:/3DTPE/analysis/TeSLA/TA 1000 unit samples');
[folderPath, folderName] = fileparts(fileparts(filePath));
load(strcat(filePath, fileName));
laneNum = 8;
totalBands = 0;
for i = 1:numel(bandStat)
    totalBands = totalBands + bandStat(i).count;
end
    
bandSizeInfo = zeros(totalBands, 3);
multiCount = 0;
for p = 1:numel(bandStat)
    bandSizeInfo(p, 3) = bandStat(p).bandPos(1);
    individualSize = bandStat(p).bandSize;
    if individualSize == 0
        individualSize = 0.00001;
    end
    bandSizeInfo(p, 2) = individualSize;
    if bandStat(p).count > 1
        for q = 2:bandStat(p).count
            multiCount = multiCount + 1;
            bandSizeInfo(numel(bandStat) + multiCount, 3) = bandStat(p).bandPos(1);            
            bandSizeInfo(numel(bandStat) + multiCount, 2) = individualSize;
        end
    end
end

bandSizeInfo(:, 1) = kmeans(bandSizeInfo(:, 3), laneNum);
bandSizeInfo = sortrows(bandSizeInfo);
indexNum = bandSizeInfo(:, 1);
bandHist = hist(indexNum, unique(indexNum));
bandSize = NaN(max(bandHist), laneNum);

count = 0;
for lane = 1:laneNum
    bandNum = bandHist(lane);
    for band = 1:bandNum
        count = count + 1;
        bandSize(band, lane) = bandSizeInfo(count, 2);
    end
end

if ~exist(strcat(folderPath, '/bandSizeStat'), 'dir')
    mkdir(strcat(folderPath, '/bandSizeStat'));
end
save(strcat(folderPath, '/bandSizeStat/', folderName, ' bandSizeStat'), 'bandSize')

end

