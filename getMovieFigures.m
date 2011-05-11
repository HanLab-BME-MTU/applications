function movieData = getMovieFigures(movieData, bandWidth, varargin)

% BEGIN
movieData.figures.status = 0;

% Parse input parameters
checkMovieData = @(movieData) ...
    checkMovieDistanceTransform(movieData) && ...
    checkMoviePairTracks(movieData);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', checkMovieData);
ip.addRequired('bandWidth', @isscalar);
ip.addParamValue('minActinLifetime', 3, @isscalar);
ip.addParamValue('minSegmentsPerBin', 10, @isscalar);
ip.addParamValue('batchMode', true, @islogical);

ip.parse(movieData, bandWidth, varargin{:});
minActinLifetime = ip.Results.minActinLifetime;
minSegmentsPerBin = ip.Results.minSegmentsPerBin;
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
[trackInfos, xMap, yMap] = mpm2trackInfos(MPM, distToEdge, [0 bandWidth],...
    minActinLifetime, true);

row = trackInfos{1}(:,1);

tFirst = trackInfos{1}(:,2);
tLast = trackInfos{1}(:,3);
lifetime = tLast - tFirst + 1;

pFirst = sub2ind(size(xMap), row, tFirst);
pLast = sub2ind(size(xMap), row, tLast);

speeds = sqrt((xMap(pFirst) - xMap(pLast)).^2 + ...
    ((yMap(pFirst) - yMap(pLast)).^2)) ./ (lifetime - 1);
speeds = speeds * pixelSize * 60 / timeInterval;

%% Find what is the average Actin speed in the vicinity of each adhesion,

tmp = vertcat(segments{:}); %#ok<USENS>
tmp = num2cell(tmp(:,1:5),1);
[x1,y1,x2,y2,L] = tmp{:};
maxAdhesionLength = sqrt(max((x2 - x1).^2 + (y2 - y1).^2)) * pixelSize;
nAdhesions = max(L);
rangeLength = [0, 500:250:1000, 1500:500:(maxAdhesionLength + 500)];
nBins = numel(rangeLength) - 1;

accCount = zeros(nAdhesions, nBins);
accSpeed = zeros(nAdhesions, nBins);

for iFrame = 1:nFrames
    % Find which Actin track leaves in iFrame
    isInFrame = iFrame >= tFirst & iFrame <= tLast;
    
    if ~any(isInFrame)
        continue;
    end

    X = xMap(row(isInFrame), iFrame);
    Y = yMap(row(isInFrame), iFrame);
    ind = sub2ind(imSize, Y, X);
    speedMap = zeros(imSize);
    speedMap(ind) = speeds(isInFrame);
    
    % Find which adhesion belongs to the frame (idx)
    segmentParams = num2cell(segments{iFrame},1);
    [x1,y1,x2,y2,idx] = segmentParams{1:5};
    
    % Calculate the Actin in the vicinity of each adhesion within iFrame
    x = .5 * (x1 + x2);
    y = .5 * (y1 + y2);
    l = sqrt((x2 - x1).^2 + (y2 - y1).^2) * pixelSize;
    s = repmat(sigmaPSF * kSigma, numel(x), 1);
    t = atan2(y2 - y1, x2 - x1);

    [xRange, yRange, nzIdx] = arrayfun(@(x,y,l,s,t) ...
        segment2DSupport(x,y,l,s,t,kSigma,fliplr(imSize)), x, y, l, s, t, ...
        'UniformOutput', false);
    
    speedMapCrops = cellfun(@(xRange,yRange) speedMap(yRange,xRange), ...
        xRange, yRange, 'UniformOutput', false);
    
    % Find out which bin each adhesion falls into
    for iBin = 1:nBins
        isInBin = rangeLength(iBin) <= l & l < rangeLength(iBin+1);
        ind = idx(isInBin);
        
        accCount(ind,iBin) = accCount(ind,iBin) + 1;
        accSpeed(ind,iBin) = accSpeed(ind,iBin) + cellfun(@(crop,nzIdx) ...
            mean(nonzeros(crop(nzIdx))), speedMapCrops(isInBin), ...
            nzIdx(isInBin));
    end
end

isValid = accCount ~= 0;
accSpeed(isValid) = accSpeed(isValid) ./ accCount(isValid);

data = num2cell(accSpeed,1);
data = cellfun(@(c) nonzeros(c(~isnan(c))), data, 'UniformOutput', false);

prm = zeros(5,nBins);

for iBin = 1:nBins
    D = sort(data{iBin});
     
    if numel(D) >= minSegmentsPerBin
        prm(1,iBin) = D(floor(numel(D)/2)+1);
        prm(2,iBin) = D(floor(numel(D)/4)+1);
        prm(3,iBin) = D(floor(3 * numel(D)/4)+1);
        prm(4,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
        prm(5,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
    else
        prm(:,iBin) = NaN;
    end
     
    data{iBin} = D;
end
 
hFig = figure('Visible', 'off');
hold on;

xlabels = arrayfun(@(iBin) [num2str(rangeLength(iBin)) '-' num2str(rangeLength(iBin+1))], ...
    2:nBins, 'UniformOutput', false);

xlabels = ['Diffraction-limited', xlabels];

boxplot2({prm},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', ...
    'Actin Speed (nm/min)');

hAxes = get(hFig, 'CurrentAxes');
XTicks = get(hAxes, 'XTick');
YTicks = cellfun(@numel, data);
hold on;
[AX,~,H2] = plotyy(NaN,NaN,XTicks,YTicks);
set(AX(2),'XTickLabel',[]);
set(H2,'Color', 'r');
set(AX(2),'YColor', 'r');
set(AX(1),'YColor', 'k');
ylabel(AX(2), 'Adhesion Count', 'Color', 'r');
ylabel(AX(1), 'Actin Speed (nm/min)', 'Color', 'k');

% Perform ks-test2 between the 2 first categories
h = kstest2(data{1}, data{2});

if h
    y = max(prm(3,1:2) + prm(4,1:2));
    line(XTicks(1:2), repmat(y + 40, 1, 2) + 10, 'Color', 'k', 'LineWidth', 4);
    text(mean(XTicks(1:2)), y + 65, '*', 'FontName', 'Helvetica', 'FontSize', 24);
end

%% Saving
fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_fig4C.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%% END
movieData.figures.status = 1;
