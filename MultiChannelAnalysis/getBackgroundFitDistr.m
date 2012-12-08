%Estimates the amplitude distribution for fits at random locations in the cell background

% Francois Aguet, 11/2012

function A = getBackgroundFitDistr(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('np', 10000);
ip.addParamValue('ccpMask', 'on', @(x) any(strcmpi(x, {'on','off'})));
ip.parse(data, varargin{:});

% load cell mask
cellmask = logical(getCellMask(data));
[ny,nx] = size(cellmask);

% Determine master/slave channels
nCh = length(data.channels); % number of channels
mCh = strcmp(data.source, data.channels);

sigma = arrayfun(@(k) getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{k}), 1:nCh);
% used for boundary and initializations only, take largest sigma
w = ceil(4*max(sigma));


frameIdx = round(linspace(1, data.movieLength, 12));
nf = numel(frameIdx);

A = cell(1,nf);

parfor i = 1:nf
    k = frameIdx(i);
    
    %-----------------------------------------------
    % Generate masks
    %-----------------------------------------------
    mask = cellmask;
    if strcmpi(ip.Results.ccpMask, 'on')
        ccpMask = double(imread(data.maskPaths{k}));
        ccpMask(ccpMask~=0) = 1;
        mask = cellmask - imdilate(ccpMask, strel('disk', 1*ceil(4*sigma(mCh))));
        % mask of 'endocytically active' zone: dilation of CCP detections
        %mask = cellmask & imdilate(ccpMask, strel('disk', 1*ceil(4*sigma(mCh))));
    end
    
    % generate candidate points
    N = ip.Results.np;
    x = [];
    y = [];
    linIdx = [];
    while numel(x)<N
        xcand = (nx-2*w-1)*rand(1,N)+w+1;
        ycand = (ny-2*w-1)*rand(1,N)+w+1;
        xi = round(xcand);
        yi = round(ycand);
    
        % remove points outside of mask or within border
        linIdxCand = sub2ind([ny nx], yi, xi);
        %rmIdx = cellmask(linIdxCand)==0 | xi<=w | yi<=w | xi>nx-w | yi>ny-w;
        validIdx = mask(linIdxCand)==1 & xi>w & yi>w & xi<=nx-w & yi<=ny-w;
        x = [x xcand(validIdx)];
        y = [y ycand(validIdx)];
        linIdx = [linIdx linIdxCand(validIdx)];
    end
    x = x(1:N);
    y = y(1:N);
    linIdx = linIdx(1:N);
    
    % estimate amplitude in each channel
    A{i} = zeros(nCh,N);
    for c = 1:nCh
        frame = double(imread(data.framePaths{c}{k}));
        
        % get local min & max for initial c and A
        ww = 2*w+1;
        maxF = ordfilt2(frame, ww^2, true(ww));
        minF = ordfilt2(frame, 1, true(ww));
        pStruct = fitGaussians2D(frame, x, y, maxF(linIdx), sigma, minF(linIdx), 'Ac');
        A{i}(c,:) = pStruct.A;
    end
end
