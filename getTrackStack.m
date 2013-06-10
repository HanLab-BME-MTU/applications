% stack = getTrackStack(data, track)
%     Returns the image data around a track in a 4-sigma window as a cell array
%     containing all channels.

% Francois Aguet, Jan 26 2011 (last modified 02/03/2012)

function [stack, xa, ya] = getTrackStack(data, track, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('track', @isstruct);
ip.addParamValue('Reference', 'track', @(x) strcmpi(x, 'track') | strcmpi(x, 'frame'));
ip.addParamValue('WindowWidth', 5, @isscalar);
ip.parse(data, track, varargin{:});

nc = length(data.channels);
mCh = strcmp(data(1).source, data(1).channels);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{mCh});
w = ceil(ip.Results.WindowWidth*sigma);

% coordinate matrices
xv = track.x;
yv = track.y;

% start and end buffer sizes
if ~isempty(track.startBuffer)
    sb = numel(track.startBuffer.t);
    xv = [track.startBuffer.x xv];
    yv = [track.startBuffer.y yv];
else
    sb = 0;
end
if ~isempty(track.endBuffer)
    eb = numel(track.endBuffer.t);
    xv = [xv track.endBuffer.x];
    yv = [yv track.endBuffer.y];
else
    eb = 0;
end

% frame index
fi = track.start-sb:track.end+eb;
nf = length(fi);

stack = cell(nc,nf);

if track.nSeg==1 && strcmpi(ip.Results.Reference, 'track') % align frames to track
    xi = round(xv(mCh,:));
    yi = round(yv(mCh,:));
    % ensure that window falls within frame bounds
    x0 = xi - min([xi-1 w]);
    x1 = xi + min([data.imagesize(2)-xi w]);
    y0 = yi - min([yi-1 w]);
    y1 = yi + min([data.imagesize(1)-yi w]);
    % axes for each frame
    xa = arrayfun(@(i) x0(i):x1(i), 1:nf, 'unif', 0);
    ya = arrayfun(@(i) y0(i):y1(i), 1:nf, 'unif', 0);
else
    % window around track mean
    mu_x = round(nanmean(xv,2));
    mu_y = round(nanmean(yv,2));
    x0 = max(1, min(mu_x)-w);
    x1 = min(data.imagesize(2), max(mu_x)+w);
    y0 = max(1, min(mu_y)-w);
    y1 = min(data.imagesize(1), max(mu_y)+w);
    xa = repmat({x0:x1}, [nf 1]);
    ya = repmat({y0:y1}, [nf 1]);
end

% load all visible frames of this track and store
for c = 1:nc
    for k = 1:nf
        frame = imread(data.framePaths{c}{fi(k)});
        stack{c,k} = frame(ya{k}, xa{k});
    end
end
