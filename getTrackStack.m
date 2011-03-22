% stack = getTrackStack(data, track)
%     Returns the image data around a track in a 4-sigma window as a cell array
%     containing all channels.

% Francois Aguet, Jan 26 2011

function stack = getTrackStack(data, track)

nChannels = length(data.channels);

masterChannel = find(strcmp(data(1).source, data(1).channels));
slaveChannels = setdiff(1:nChannels, masterChannel);


sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, name2wavelength(data.markers{masterChannel}));
w = ceil(4*sigma);

% buffer with 5 frames before and after
buffer = 5;
bStart = track.start - max(1, track.start-buffer);
bEnd = min(data.movieLength, track.end+buffer) - track.end;
idx = track.start-bStart:track.end+bEnd;
nf = length(idx);

stack = cell(nChannels,nf);
    
xi = round(track.x);
yi = round(track.y);
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

