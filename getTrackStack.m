% stack = getTrackStack(data, track)
%     Returns the image data around a track in a 4-sigma window as a cell array
%     containing all channels.

% Francois Aguet, Jan 26 2011

function stack = getTrackStack(data, track, varargin)

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
masterChannel = find(strcmp(data(1).source, data(1).channels));
slaveChannels = setdiff(1:nc, masterChannel);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, name2wavelength(data.markers{masterChannel}));
w = ceil(ip.Results.WindowWidth*sigma);

% buffer with 5 frames before and after
bStart = track.start - max(1, track.start-buffer);
bEnd = min(data.movieLength, track.end+buffer) - track.end;
idx = track.start-bStart:track.end+bEnd;
nf = length(idx);

stack = cell(nc,nf);

if strcmpi(ip.Results.Reference, 'track')
    xi = round(track.x);
    yi = round(track.y);
else
    xi = ones(1,nf)*round(mean([track.x]));
    yi = ones(1,nf)*round(mean([track.y]));
end

xi = [xi(1)*ones(1,bStart) xi xi(end)*ones(1,bEnd)];
yi = [yi(1)*ones(1,bStart) yi yi(end)*ones(1,bEnd)];
% load all visible frames of this track and store
for c = [masterChannel slaveChannels]
    tifFiles = dir([data.channels{c} '*.tif*']);
    tifFiles = tifFiles(idx);
    for k = 1:nf
        frame = imread([data.channels{c} tifFiles(k).name]);
        stack{c,k} = frame(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    end
end

