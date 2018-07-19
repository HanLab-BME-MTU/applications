function [ratios proportions] = computeRatioIntensity(nChannels, bandWidth, normalized)

if nChannels <= 1
    error('nChannels needs to be greater than 1.');
end

% Ask for mask
[fileName,pathName] = uigetfile({'*.tif'},'Select mask file');

if ~ischar(fileName) || ~ischar(pathName)
    disp('Abord');
    return;
end

mask = imread(fullfile(pathName, fileName));

fileNames = cell(nChannels,1);
pathNames = cell(nChannels,1);

for iChannel = 1:nChannels
    [fileName,pathName] = uigetfile({'*.tif'},['Select image file for channel ' num2str(iChannel)]);
    
    if ~ischar(fileName) || ~ischar(pathName)
        disp('Abord');
        return;
    end

    fileNames{iChannel} = fileName;
    pathNames{iChannel} = pathName;
end

dist = bwdist(1 - double(mask));

maxDist = max(dist(:));

bands = 0:bandWidth:maxDist;

% Compute ratios

ratios = zeros(numel(bands) - 1, nChannels * (nChannels - 1) / 2);

iRatio = 1;

for iChannel = 1:nChannels-1    
    I = imread(fullfile(pathNames{iChannel}, fileNames{iChannel}));
    for jChannel = iChannel+1:nChannels
        J = imread(fullfile(pathNames{jChannel}, fileNames{jChannel}));
        
        for iBand = 1:numel(bands)-1
            ratios(iBand, iRatio) = ...
                mean(I(dist > bands(iBand) & dist <= bands(iBand+1))) / ...
                mean(J(dist > bands(iBand) & dist <= bands(iBand+1)));
        end
        
        iRatio = iRatio + 1;
    end
end

% Compute proportions

proportions = zeros(numel(bands) - 1, nChannels);

for iChannel = 1:nChannels
    I = imread(fullfile(pathNames{iChannel}, fileNames{iChannel}));

    for iBand = 1:numel(bands)-1
        proportions(iBand, iChannel) = mean(I(dist > bands(iBand) & dist <= bands(iBand+1)));
    end
end

figure, bar(ratios);
xlabel('Bands (pixel)');
ylabel('Average intensity');
title('Ratios');

if normalized
    proportions = bsxfun(@rdivide,proportions,sum(proportions,2));
end

figure, bar(proportions,'stack');
xlabel('Bands (pixel)');
ylabel('Average intensity');
title('Proportions');