function computeRatioIntensity(nChannels, bandWidth, normalized)

% Ask for mask
[fileName,pathName] = uigetfile({'*.tif'},'Select mask file');

mask = imread(fullfile(pathName, fileName));

fileNames = cell(nChannels,1);
pathNames = cell(nChannels,1);

for iChannel = 1:nChannels
    [fileName,pathName] = uigetfile({'*.tif'},['Select image file for channel ' num2str(iChannel)]);
    
    fileNames{iChannel} = fileName;
    pathNames{iChannel} = pathName;
end

dist = bwdist(1 - double(mask));

maxDist = max(dist(:));

bands = 0:bandWidth:maxDist;

ratios = zeros(numel(bands) - 1, nChannels);

for iChannel = 1:nChannels
    I = imread(fullfile(pathNames{iChannel}, fileNames{iChannel}));

    for iBand = 1:numel(bands)-1
        ratios(iBand, iChannel) = mean(I(dist > bands(iBand) & dist <= bands(iBand+1)));
    end
end

if normalized
    ratios = bsxfun(@rdivide,ratios,sum(ratios,2));
end

bar(ratios,'stack');
xlabel('Bands (pixel)');
ylabel('Average intensity');