function movieData = getMovieFigures(movieData,batchMode)

assert(checkMovieDistanceTransform(movieData));
assert(checkMoviePairTracks(movieData));

movieData.figures.status = 0;

movieData.figures.directory = fullfile(movieData.analysisDirectory, 'figures');

% Create output directory
if ~exist(movieData.figures.directory, 'dir')
    mkdir(movieData.figures.directory);
end

nFrames = movieData.nImages(1);
imSize = movieData.imSize;
pixelSize = movieData.pixelSize_nm;
timeInterval = movieData.timeInterval_s;
sigmaPSF = movieData.particleDetection.params.sigmaPSF;
kSigma = movieData.particleDetection.params.kSigma;

%% Load data

load(fullfile(movieData.pairTracks.directory, ['CCParams_iter=' ...
    num2str(movieData.pairTracks.params.maxIter-1) '_.mat']));

% Check that Actin tracking has been done
load(fullfile(movieData.imageDirectory, movieData.channelDirectory{2}, ...
    '..', 'analysis', 'tack', 'mpm.mat'));

% Read the list of distance transforms
distToEdgePath = movieData.distanceTransform.directory;
distToEdgeFiles = dir([distToEdgePath filesep '*.mat']);
    
% Read distance transforms
distToEdge = zeros(movieData.imSize(1),movieData.imSize(2),nFrames);
    
for iFrame = 1:nFrames
    fileName = fullfile(distToEdgePath, distToEdgeFiles(iFrame).name);
    tmp = load(fileName);
    distToEdge(:,:,iFrame) = tmp.distToEdge * pixelSize;
end

%% Rearrange Actin tracks and get Actin speed
[trackInfos, xMap, yMap] = mpm2trackInfos(MPM,distToEdge,[0 +Inf],kSigma,...
    fliplr(imSize));

row = trackInfos{1}(:,1);

tFirst = trackInfos{1}(:,2);
tLast = trackInfos{1}(:,3);
lifetime = tLast - tFirst + 1;

pFirst = sub2ind(size(xMap), row, tFirst);
pLast = sub2ind(size(xMap), row, tLast);

speeds = sqrt((xMap(pFirst) - xMap(pLast)).^2 + ...
    ((yMap(pFirst) - yMap(pLast)).^2)) ./ (lifetime - 1);
speeds = speeds * pixelSize * 60 / timeInterval;

%% Compute the correlation
nSegmentsPerFrame = cellfun(@(c) size(c,1), segments); %#ok<USENS>

ppLast = cumsum(nSegmentsPerFrame);
ppFirst = ppLast - nSegmentsPerFrame + 1;

actinSpeedPerSegment = nan(sum(nSegmentsPerFrame), 1);
lengthPerSegment = nan(sum(nSegmentsPerFrame), 1);

for iFrame = min(tFirst):max(tLast)
    % Find which Actin track leaves in iFrame
    isInFrame = iFrame >= tFirst & iFrame <= tLast;
    
    % Find the location of each track in the frame
    X = xMap(row(isInFrame),iFrame);
    Y = yMap(row(isInFrame),iFrame);
    ind = sub2ind(imSize,Y,X);
    speedMap = zeros(imSize);
    speedMap(ind) = speeds(isInFrame);
    
    segmentParams = num2cell(segments{iFrame},1);
    [x1,y1,x2,y2] = segmentParams{:};
    x = .5 * (x1 + x2);
    y = .5 * (y1 + y2);
    l = sqrt((x2 - x1).^2 + (y2 - y1).^2);
    s = repmat(sigmaPSF * kSigma, numel(x), 1);
    t = atan2(y2 - y1, x2 - x1);
    
    [xRange, yRange, nzIdx] = arrayfun(@(x,y,l,s,t) ...
        segment2DSupport(x,y,l,s,t,kSigma,fliplr(imSize)), x, y, l, s, t, ...
        'UniformOutput', false);
    
    speedMapCrops = cellfun(@(xRange,yRange) speedMap(yRange,xRange), ...
        xRange, yRange,  'UniformOutput', false);
    
    actinSpeedPerSegment(ppFirst(iFrame):ppLast(iFrame)) = ...
        cellfun(@(crop,nzIdx) mean(nonzeros(crop(nzIdx))), ...
        speedMapCrops, nzIdx);
    
    lengthPerSegment(ppFirst(iFrame):ppLast(iFrame)) = l * pixelSize;
end

%% Output
isValid = ~isnan(actinSpeedPerSegment);
lengthPerSegment = lengthPerSegment(isValid);
actinSpeedPerSegment = actinSpeedPerSegment(isValid);

range = [0:250:1000, 1500:500:max(lengthPerSegment)];

nBins = numel(range)-1;

prm = zeros(5,nBins);

for iBin = 1:nBins
    isInBin = lengthPerSegment >= range(iBin) & lengthPerSegment < range(iBin+1);
    
    data = sort(actinSpeedPerSegment(isInBin));
    
    if numel(data)
        prm(1,iBin) = data(floor(numel(data)/2)+1);
        prm(2,iBin) = data(floor(numel(data)/4)+1);
        prm(3,iBin) = data(floor(3 * numel(data)/4)+1);
        prm(4,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
        prm(5,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
    else
        prm(:,iBin) = NaN;
    end
end

hFig = figure('Visible', 'off');
hold on;

xlabels = arrayfun(@(iBin) [num2str(range(iBin)) '-' num2str(range(iBin+1))], ...
    1:nBins, 'UniformOutput', false);

boxplot2({prm},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', 'Actin Speed (nm/min)');

fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_fig4C.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%% END
movieData.figures.status = 1;
