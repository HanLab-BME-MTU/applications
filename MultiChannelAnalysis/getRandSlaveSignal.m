function track = getRandSlaveSignal(data, track, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('track', @isstruct);
ip.parse(data, track, varargin{:});

% load cell mask
cellmask = logical(getCellMask(data));
[ny,nx] = size(cellmask);

% Determine master/slave channels
nCh = length(data.channels); % number of channels
mCh = strcmp(data.source, data.channels);
sCh = setdiff(1:nCh, mCh);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{sCh});
w = ceil(4*sigma);

b = 5;
frameIdx = track.f(1)-b:track.f(end)+b;

nf = numel(frameIdx);
A = zeros(1,nf);
A_pstd = zeros(1,nf);
c = zeros(1,nf);
c_pstd = zeros(1,nf);

for i = 1:nf;
    k = frameIdx(i);

    frame = double(imread(data.framePaths{2}{k}));
    
    pointFound = false;
    while ~pointFound
    
        % generate candidate points
        x = (nx-2*w-1)*rand+w+1;
        y = (ny-2*w-1)*rand+w+1;
        xi = round(x);
        yi = round(y);
    
        % remove points outside of mask or within border
        linIdx = sub2ind([ny nx], yi, xi);
        pointFound = cellmask(linIdx)==1 & xi>w & yi>w & xi<=nx-w & yi<=ny-w;
    end
    
    % get local min & max for initial c and A
    window = frame(yi-w:yi+w, xi-w:xi+w);

    maxF = max(window(:));
    minF = min(window(:));
    
    pstruct = fitGaussians2D(frame, x, y, maxF, sigma, minF, 'Ac');
    A(i) = pstruct.A;
    A_pstd(i) = pstruct.A_pstd;
    c(i) = pstruct.c;
    c_pstd(i) = pstruct.c_pstd;
end
track.A(2,:) = A(b+1:end-b);
track.A_pstd(2,:) = A_pstd(b+1:end-b);
track.c(2,:) = c(b+1:end-b);
track.c_pstd(2,:) = c_pstd(b+1:end-b);

track.startBuffer.A(2,:) = A(1:b);
track.startBuffer.A_pstd(2,:) = A_pstd(1:b);
track.startBuffer.c(2,:) = c(1:b);
track.startBuffer.c_pstd(2,:) = c_pstd(1:b);

track.endBuffer.A(2,:) = A(end-b+1:end);
track.endBuffer.A_pstd(2,:) = A_pstd(end-b+1:end);
track.endBuffer.c(2,:) = c(end-b+1:end);
track.endBuffer.c_pstd(2,:) = c_pstd(end-b+1:end);
