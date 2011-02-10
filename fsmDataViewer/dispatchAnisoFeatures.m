function dataLayer = dispatchAnisoFeatures(fileList, nFrames)

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
    t = featuresInfo(iFrame).theta(:,1);
    ct = cos(t);
    st = sin(t);

    dataLayer{iFrame} = [y, x, 5 * st, 5 * ct];
end
