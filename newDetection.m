function newDetection(data, frameIdx)
if nargin<2
    frameIdx = 1;
end

% master channel
mCh = strcmp(data.channels, data.source);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{mCh});

% xa = (0:data.imagesize(2)-1)*data.pixelSize/data.M;
% ya = (0:data.imagesize(1)-1)*data.pixelSize/data.M;

for k = frameIdx%:data.movieLength
    img = double(imread(data.framePaths{mCh}{k}));
    %img = img(265:300,880:980);
    %[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma);
    
    % W: W1, W2, W3, W4, A1 %%%%%
    %W = awt(img, 4);
    %[wtStruct, wtMask] = detectSpotsWT(img);
    wtMask = double(imread(data.maskPaths{k}));
    [pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma);
    
    % intensity for clusters
    %[labels nComp] = bwlabel(mask);   
end

rgbMask = rgbOverlay(img, {mask, wtMask}, {[0 1 0], [1 0 0]});
% rgbMask = rgbOverlay(img, {logMask, mask}, {[0 1 0],[1 0 0]});

figure('Units', 'Normalized', 'Position', [0 0 1 1], 'PaperPositionMode', 'auto');
ha1 = axes('Position', [0 0.55 1 0.45]);
imagesc(rgbMask); axis image;
hold on;
% initial position/local maxima
%[lmy, lmx] = find(imgLM~=0);
%plot(lmx,lmy, 'wx');

xr = [pstruct.x];
yr = [pstruct.y];
% plot(xr, yr, 'ro');

plot(xr(pstruct.isPSF), yr(pstruct.isPSF), 's', 'Color', 0.6*[1 1 1], 'MarkerSize', 10);
plot(xr(~pstruct.isPSF), yr(~pstruct.isPSF), 'rx', 'MarkerSize', 10);


ha2 = axes('Position', [0 0.05 1 0.45]);
imagesc(img); axis image; colormap(gray(256));

linkaxes([ha1 ha2]);

