function dataLayer = dispatchSegment2D(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

load(fileList{1});

if ~exist('segments', 'var')
    error('Unable to find segments.');
end

if numel(segments) ~= nFrames
    error('Size of segments does not match the number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    params = num2cell(segments{iFrame},1);
    [X, Y, ~, L, S, T, ~] = params{:};
    L = L/2;
    CT = cos(T);
    ST = sin(T);
    
    L1 = repmat([X,Y],1,2);
    L2 = repmat([X,Y],1,2);
    
    L1(:,1) = L1(:,1) + L .* CT; L1(:,2) = L1(:,2) + L .* ST;
    L1(:,3) = L1(:,3) - L .* CT; L1(:,4) = L1(:,4) - L .* ST;
    L2(:,1) = L2(:,1) - S .* ST; L2(:,2) = L2(:,2) + S .* CT;
    L2(:,3) = L2(:,3) + S .* ST; L2(:,4) = L2(:,4) - S .* CT;
    
    dataLayer{iFrame} = vertcat(L1,L2);
end