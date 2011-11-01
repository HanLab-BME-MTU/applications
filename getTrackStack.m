% stack = getTrackStack(data, track)
%     Returns the image data around a track in a 4-sigma window as a cell array
%     containing all channels.

% Francois Aguet, Jan 26 2011 (last modified 10/26/2011)

function [stack, xa, ya] = getTrackStack(data, track, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('track', @isstruct);
ip.addParamValue('Reference', 'track', @(x) strcmpi(x, 'track') | strcmpi(x, 'frame'));
ip.addParamValue('WindowWidth', 8, @isscalar);
ip.addParamValue('Buffer', 5, @isscalar);
ip.addParamValue('Channels', []);
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
    xi = round(xv{1}(mCh,:));
    yi = round(yv{1}(mCh,:));
    % ensure that window falls within frame bounds
    
    x0 = xi - min([xi-1 w]);
    x1 = xi + min([nx-xi w]);
    y0 = yi - min([yi-1 w]);
    y1 = yi + min([ny-yi w]);
    xa = arrayfun(@(i) x0(i):x1(i), 1:nf, 'UniformOutput', false);
    ya = arrayfun(@(i) y0(i):y1(i), 1:nf, 'UniformOutput', false);
else
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
        frame = imread(data.framePaths{c}{k});
        stack{c,k} = frame(ya{1}, xa{k});
    end
end
