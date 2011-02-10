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

X = [featuresInfo(:).xCoord]; X = X(:,1);
Y = [featuresInfo(:).yCoord]; Y = Y(:,1);
SX = [featuresInfo(:).stdAlong]; SX = SX(:,1);
SY = [featuresInfo(:).stdAside]; SY = SY(:,1);
T = [featuresInfo(:).xCoord]; T = T(:,1);
CT = cos(T);
ST = sin(T);

L1 = repmat([X,Y],1,2);
L2 = repmat([X,Y],1,2);

L1(:,1) = L1(:,1) + SX .* CT; L1(:,2) = L1(:,2) - SX .* ST;
L1(:,3) = L1(:,3) - SX .* CT; L1(:,4) = L1(:,4) + SX .* ST;
L2(:,1) = L2(:,1) - SY .* ST; L2(:,2) = L2(:,2) - SY .* CT;
L2(:,3) = L2(:,3) + SY .* ST; L2(:,4) = L2(:,4) + SY .* CT;

dataLayer{iFrame} = vertcat(L1,L2);
