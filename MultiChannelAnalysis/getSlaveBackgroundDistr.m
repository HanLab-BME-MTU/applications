function A = getSlaveBackgroundDistr(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('np', 10000);
ip.parse(data, varargin{:});


% load cell mask
cellmask = logical(getCellMask(data));
[ny,nx] = size(cellmask);

% Determine master/slave channels
nCh = length(data.channels); % number of channels
mCh = strcmp(data.source, data.channels);
sCh = setdiff(1:nCh, mCh);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{sCh});
w = ceil(4*sigma);


frameIdx = round(linspace(1, data.movieLength, 12));
nf = numel(frameIdx);

A = cell(1,nf);
parfor i = 1:nf;
    k = frameIdx(i);
    
    %-----------------------------------------------
    % Generate masks
    %-----------------------------------------------
    % load CCP mask and dilate
    ccpMask = double(imread(data.maskPaths{k}));
    ccpMask(ccpMask~=0) = 1;
    ccpMask = imdilate(ccpMask, strel('disk', 1*w));
    % mask = maskInt-imdilate(mask, strel('disk', w));
    ccpMask = ccpMask & cellmask; % endocytically active zone (EAZ)
    
    frame = double(imread(data.framePaths{2}{k}));
        
    % generate candidate points
    x = (nx-2*w-1)*rand(1,ip.Results.np)+w+1;
    y = (ny-2*w-1)*rand(1,ip.Results.np)+w+1;
    xi = round(x);
    yi = round(y);
    
    % remove points outside of mask or within border
    linIdx = sub2ind([ny nx], yi, xi);
    rmIdx = cellmask(linIdx)==0 | xi<=w | yi<=w | xi>nx-w | yi>ny-w;
    x(rmIdx) = [];
    y(rmIdx) = [];
    linIdx(rmIdx) = [];
    
    % get local min & max for initial c and A
    ww = 2*w+1;
    maxF = ordfilt2(frame, ww^2, true(ww));
    minF = ordfilt2(frame, 1, true(ww));
    
    pStruct = fitGaussians2D(frame, x, y, maxF(linIdx), sigma, minF(linIdx), 'Ac');
    A{i} = pStruct.A;
end