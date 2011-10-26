% stack = getTrackStack(data, track)
%     Returns the image data around a track in a 4-sigma window as a cell array
%     containing all channels.

% Francois Aguet, Jan 26 2011 (last modified 10/26/2011)

function [stack, dx, dy] = getTrackStack(data, track, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('track', @isstruct);
ip.addParamValue('Reference', 'track', @(x) strcmpi(x, 'track') | strcmpi(x, 'frame'));
ip.addParamValue('WindowWidth', 4, @isscalar);
ip.addParamValue('Buffer', 5, @isscalar);
ip.parse(data, track, varargin{:});
buffer = ip.Results.Buffer;

nc = length(data.channels);
mCh = strcmp(data(1).source, data(1).channels);
ny = data.imagesize(1);
nx = data.imagesize(2);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, name2wavelength(data.markers{mCh}));
w = ceil(ip.Results.WindowWidth*sigma);

% buffer with 5 frames before and after
bStart = track.start - max(1, track.start-buffer);
bEnd = min(data.movieLength, track.end+buffer) - track.end;
idx = track.start-bStart:track.end+bEnd;
nf = length(idx);

stack = cell(nc,nf);

ns = numel(track.x);
xv = arrayfun(@(s) [track.startBuffer.x{s} track.x{s} track.endBuffer.x{s}], 1:ns, 'UniformOutput', false);
yv = arrayfun(@(s) [track.startBuffer.y{s} track.y{s} track.endBuffer.y{s}], 1:ns, 'UniformOutput', false);

mu_x = cellfun(@(s) round(mean(s,2)), xv, 'UniformOutput', false);
mu_x = [mu_x{:}];
mu_x = mu_x(mCh,:);
mu_y = cellfun(@(s) round(mean(s,2)), yv, 'UniformOutput', false);
mu_y = [mu_y{:}];
mu_y = mu_y(mCh,:);

if ns==1 && strcmpi(ip.Results.Reference, 'track') % align frames to track
    xi = round(xv{1});
    yi = round(yv{1});
    % ensure that window falls within frame bounds
    x0 = xi - min(min(xi,[],2)-1,w,[],2);
    x1 = xi + min(nx-max(xi,[],2),w,[],2);
    y0 = yi - min(min(yi,[],2)-1,w,[],2);
    y1 = yi + min(ny-max(yi,[],2),w,[],2);
    dx = xv{1}-xi;
    dy = yv{1}-yi;
else
    x0 = repmat(max(1, min(mu_x)-w), [1 nf]);
    x1 = repmat(min(data.imagesize(2), max(mu_x)+w), [1 nf]);
    y0 = repmat(max(1, min(mu_y)-w), [1 nf]);
    y1 = repmat(min(data.imagesize(1), max(mu_y)+w), [1 nf]);
    dx = [];
    dy = [];
end

% load all visible frames of this track and store
for c = 1:nc
    for k = 1:nf
        frame = imread(data.framePaths{c}{k});
        stack{c,k} = frame(y0(k):y1(k), x0(k):x1(k));
    end
end
