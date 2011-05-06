function movieData = getMovieFigures(movieData, varargin)

% BEGIN
movieData.figures.status = 0;

% Parse input parameters
checkMovieData = @(movieData) ...
    checkMovieDistanceTransform(movieData) && ...
    checkMoviePairTracks(movieData);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', checkMovieData);
ip.addParamValue('batchMode', true, @islogical);

ip.parse(movieData, varargin{:});
batchMode = ip.Results.batchMode;

movieData.figures.directory = fullfile(movieData.analysisDirectory, 'figures');

% Create output directory
if ~exist(movieData.figures.directory, 'dir')
    mkdir(movieData.figures.directory);
end

nFrames = movieData.nImages(1);
imSize = movieData.imSize;
pixelSize = movieData.pixelSize_nm;
timeInterval = movieData.timeInterval_s;
sigmaPSF = movieData.particleDetection.params.required.sigmaPSF;
kSigma = movieData.particleDetection.params.optional.kSigma;

%% Load data

load(fullfile(movieData.pairTracks.directory, 'ClassifiedSegments.mat'));

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
distancePerSegment = nan(sum(nSegmentsPerFrame), 1);

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
    [x1,y1,x2,y2] = segmentParams{1:4};
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
    
    % Actin Speed
    actinSpeedPerSegment(ppFirst(iFrame):ppLast(iFrame)) = ...
        cellfun(@(crop,nzIdx) mean(nonzeros(crop(nzIdx))), ...
        speedMapCrops, nzIdx);
    
    % Length
    lengthPerSegment(ppFirst(iFrame):ppLast(iFrame)) = l * pixelSize;

    % Distance
    xi = round(x);
    yi = round(y);
    isInside = xi > 0 & yi > 0 & xi <= imSize(2) & yi <= imSize(1);
    ind = sub2ind(imSize, yi(isInside), xi(isInside));
    ind = ind + prod(imSize) * (iFrame - 1);
    dist = nan(size(xi,1),1);
    dist(isInside) = distToEdge(ind);    
    distancePerSegment(ppFirst(iFrame):ppLast(iFrame)) = dist;
end

%% Output Length / Speed
isValid = ~(isnan(actinSpeedPerSegment) | isnan(distancePerSegment));
lengthPerSegment = lengthPerSegment(isValid);
actinSpeedPerSegment = actinSpeedPerSegment(isValid);
distancePerSegment = distancePerSegment(isValid);

range = [0, 500:250:1000, 1500:500:max(lengthPerSegment)];

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
    2:nBins, 'UniformOutput', false);

xlabels = ['Diffraction-limited', xlabels];

boxplot2({prm},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', ...
    'Actin Speed (nm/min)');

fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_fig4C.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%% Output Distance / Length
range = [0:100:2000, 3000:1000:max(distancePerSegment)];
nBins = numel(range)-1;

prm = zeros(5, nBins);

for iBin = 1:nBins
    isInBin = distancePerSegment >= range(iBin) & distancePerSegment < range(iBin+1);
    
    data = sort(lengthPerSegment(isInBin));
    
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

boxplot2({prm},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', ...
    'Adhesion Lengh (nm)');

fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_suppFig4C.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);


%% END
movieData.figures.status = 1;
