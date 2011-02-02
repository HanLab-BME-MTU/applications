function dataLayer = dispatchKhuloudFeatures(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

load(fileList{1});

if ~exist('featuresInfo', 'var')
    error('Unable to find feature info.');
end

if numel(featuresInfo) ~= nFrames
    error('Size of feature info does not match the number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    x = featuresInfo(iFrame).xCoord(:,1);
    y = featuresInfo(iFrame).yCoord(:,1);

    dataLayer{iFrame} = [x, y];
end
