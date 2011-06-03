function movieData = getMovieFigures(movieData, bandWidth, varargin)

% BEGIN
movieData.figures.status = 0;

% Parse input parameters
checkMovieData = @(movieData) ...
    checkMovieDistanceTransform(movieData) && ...
    checkMoviePairTracks(movieData) && ...
    checkMovieProtrusionSamples(movieData) && ...
    checkMovieLabels(movieData);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', checkMovieData);
ip.addRequired('bandWidth', @isscalar);
ip.addParamValue('minActinLifetime', 3, @isscalar);
ip.addParamValue('minSegmentsPerBin', 10, @isscalar);
ip.addParamValue('alpha', 1e-5, @isscalar);
ip.addParamValue('batchMode', true, @islogical);

ip.parse(movieData, bandWidth, varargin{:});
minActinLifetime = ip.Results.minActinLifetime;
minSegmentsPerBin = ip.Results.minSegmentsPerBin;
alpha = ip.Results.alpha;
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

% Load adhesions
load(fullfile(movieData.pairTracks.directory, 'ClassifiedSegments.mat'));

% Load fsm MPM
load(fullfile(movieData.imageDirectory, movieData.channelDirectory{2}, ...
    '..', 'analysis', 'tack', 'mpm.mat'));

% Load distance transforms
directory = movieData.distanceTransform.directory;
files = dir([directory filesep '*.mat']);
    
D = zeros(movieData.imSize(1),movieData.imSize(2),nFrames);
    
for iFrame = 1:nFrames
    fileName = fullfile(directory, files(iFrame).name);
    tmp = load(fileName);
    D(:,:,iFrame) = tmp.distToEdge * pixelSize;
end

% Load windowing labels
directory = movieData.labels.directory;
files = dir([directory filesep '*.tif']);

L = zeros(movieData.imSize(1),movieData.imSize(2),nFrames);

for iFrame = 1:nFrames
    fileName = fullfile(directory, files(iFrame).name);
    L(:,:,iFrame) = imread(fileName);
end

% Load protrusion samples
load(fullfile(movieData.protrusionSamples.directory, ...
    movieData.protrusionSamples.fileName));

%% Rearrange Actin tracks and get Actin speed
[trackInfos, xMap, yMap] = mpm2trackInfos(MPM, D, [0 bandWidth],...
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
[x1,y1,x2,y2,l] = tmp{:};
maxAdhesionLength = sqrt(max((x2 - x1).^2 + (y2 - y1).^2)) * pixelSize;
nAdhesions = max(l);
%rangeLength = [0, 500:250:1000, 1500:500:(maxAdhesionLength + 500)];
rangeLength = [0, 500, (maxAdhesionLength + 500)];
%rangeDistance = [0:100:1000, 1500:2000:bandWidth];
rangeDistance = [0 1500 bandWidth];

nLengthBins = numel(rangeLength) - 1;
nDistanceBins = numel(rangeDistance) - 1;

accSpeedInLengthBin = zeros(nAdhesions, nLengthBins);
accSpeedInDistanceBin = zeros(nAdhesions, nDistanceBins);
accLengthInDistanceBin = zeros(nAdhesions, nDistanceBins);

accLengthBinCount = zeros(nAdhesions, nLengthBins);
accDistanceBinCount = zeros(nAdhesions, nDistanceBins);

for iFrame = 1:nFrames-1
    % Find which Actin track leaves in iFrame
    isInFrame = iFrame >= tFirst & iFrame <= tLast;
    
    if ~any(isInFrame)
        continue;
    end

    X = xMap(row(isInFrame), iFrame);
    Y = yMap(row(isInFrame), iFrame);
    ind = sub2ind(imSize, Y, X);
    speedMap = nan(imSize);
    speedMap(ind) = speeds(isInFrame);

    % Disable Actin Speed when the cell is NOT protruding
    protrusionMask = ismember(L(:,:,iFrame), find(protrusionSamples.states(:,iFrame) == 2));
    speedMap(~protrusionMask) = NaN;
    
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
    
    % Find out which length bin each adhesion falls into
    for iBin = 1:nLengthBins
        isInBin = rangeLength(iBin) <= l & l < rangeLength(iBin+1);
        ind = idx(isInBin);
        
        meanSpeedPerAdh = cellfun(@(crop,nzIdx) nanmean(crop(nzIdx)), ...
            speedMapCrops(isInBin), nzIdx(isInBin));
        
        isValid = ~isnan(meanSpeedPerAdh);
        
        accLengthBinCount(ind, iBin) = accLengthBinCount(ind, iBin) + isValid;
        accSpeedInLengthBin(ind(isValid), iBin) = accSpeedInLengthBin(ind(isValid), iBin) + ...
            reshape(meanSpeedPerAdh(isValid), numel(meanSpeedPerAdh(isValid)), 1);
    end
    
    % Find out which distance bin each adhesion falls into
    indPx = sub2ind(imSize, round(y), round(x));
    indPx = prod(imSize) * (iFrame - 1) + indPx;
    for iBin = 1:nDistanceBins
        
        isInBin = rangeDistance(iBin) <= D(indPx) & D(indPx) < rangeDistance(iBin+1);
        ind = idx(isInBin);
        
        meanSpeedPerAdh = cellfun(@(crop,nzIdx) nanmean(crop(nzIdx)), ...
            speedMapCrops(isInBin), nzIdx(isInBin));
        
        isValid = ~isnan(meanSpeedPerAdh);
        
        accDistanceBinCount(ind, iBin) = accDistanceBinCount(ind, iBin) + isValid;
        accLengthInDistanceBin(ind(isValid), iBin) = accLengthInDistanceBin(ind(isValid), iBin) + ...
            reshape(l(isValid), numel(l(isValid)), 1);
        accSpeedInDistanceBin(ind(isValid), iBin) = accSpeedInDistanceBin(ind(isValid), iBin) + ...
            reshape(meanSpeedPerAdh(isValid), numel(meanSpeedPerAdh(isValid)), 1);
    end
end

isValid = accLengthBinCount ~= 0;
accSpeedInLengthBin(isValid) = accSpeedInLengthBin(isValid) ./ ...
    accLengthBinCount(isValid);

isValid = accDistanceBinCount ~= 0;
accLengthInDistanceBin(isValid) = accLengthInDistanceBin(isValid) ./ ...
    accDistanceBinCount(isValid);
accSpeedInDistanceBin(isValid) = accSpeedInDistanceBin(isValid) ./ ...
    accDistanceBinCount(isValid);

%% Generate Actin speed in function of adhesion length

data = num2cell(accSpeedInLengthBin,1);
data = cellfun(@(c) nonzeros(c(~isnan(c))), data, 'UniformOutput', false);

fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_DATA.mat']);
save(fileName, 'data');

prm = NaN(5, nLengthBins);

bins = 1:nLengthBins;
isNotEmptyBin = cellfun(@(c) numel(c) >= minSegmentsPerBin, data);
maxNonEmptyBin = max(bins(isNotEmptyBin));

for iBin = bins(isNotEmptyBin)
    sortedData = sort(data{iBin});
     
    prm(1,iBin) = sortedData(floor(numel(sortedData)/2)+1);
    prm(2,iBin) = sortedData(floor(numel(sortedData)/4)+1);
    prm(3,iBin) = sortedData(floor(3 * numel(sortedData)/4)+1);
    prm(4,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
    prm(5,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
     
    data{iBin} = sortedData;
end
 
hFig = figure('Visible', 'off');
hold on;

xlabels = arrayfun(@(iBin) [num2str(rangeLength(iBin)) '-' num2str(rangeLength(iBin+1))], ...
    2:maxNonEmptyBin, 'UniformOutput', false);

xlabels = ['Diffraction-limited', xlabels];

boxplot2({prm(:, 1:maxNonEmptyBin)},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', ...
    'Actin Speed (nm/min)');

hAxes = get(hFig, 'CurrentAxes');
XTicks = get(hAxes, 'XTick');
yTicks = 0:100:((ceil(max(prm(3,:) + prm(4,:)) / 100) + 1) * 100);
set(hAxes,'YTickLabel', num2str(yTicks'));
set(hAxes,'YTick',yTicks');
set(hAxes,'YLim', [0 yTicks(end) + 1]);

% Perform t-test between the 2 first categories
assert(all(~isnan(data{1})) && all(~isnan(data{2})));
m1 = mean(data{1});
m2 = mean(data{2});
v1 = var(data{1});
v2 = var(data{2});
n1 = length(data{1});
n2 = length(data{2});
t = (m1 - m2) / sqrt((((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)) * (1/n1 + 1/n2));
h = (1 - tcdf(t, n1 + n2 - 2)) < alpha;

disp([getDirFromPath(movieData.imageDirectory) '   p-value = ' num2str(1 - tcdf(t, n1 + n2 - 2)) '   n1 = ' num2str(n1) '   n2 = ' num2str(n2) '   meanSpeed(bin1) = ' num2str(m1)])

if h
    y = max(prm(3,1:2) + prm(4,1:2));
    line(XTicks(1:2), repmat(y + 40, 1, 2) + 10, 'Color', 'k', 'LineWidth', 4);
    text(mean(XTicks(1:2)), y + 65, '*', 'FontName', 'Helvetica', 'FontSize', 24);
end

% Saving
fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_fig4C.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%% Generate Actin speed in function of distance to edge

data = num2cell(accSpeedInDistanceBin, 1);
data = cellfun(@(c) nonzeros(c(~isnan(c))), data, 'UniformOutput', false);

fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_ActinSpeedVSdistanceToEdge_DATA.mat']);
save(fileName, 'data');

prm = NaN(5, nDistanceBins);

bins = 1:nDistanceBins;
isNotEmptyBin = cellfun(@(c) numel(c) >= minSegmentsPerBin, data);
maxNonEmptyBin = max(bins(isNotEmptyBin));

for iBin = bins(isNotEmptyBin)
    sortedData = sort(data{iBin});
     
    prm(1,iBin) = sortedData(floor(numel(sortedData)/2)+1);
    prm(2,iBin) = sortedData(floor(numel(sortedData)/4)+1);
    prm(3,iBin) = sortedData(floor(3 * numel(sortedData)/4)+1);
    prm(4,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
    prm(5,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
     
    data{iBin} = sortedData;
end
 
hFig = figure('Visible', 'off');
hold on;

xlabels = arrayfun(@(iBin) [num2str(rangeDistance(iBin)) '-' num2str(rangeDistance(iBin+1))], ...
    1:maxNonEmptyBin, 'UniformOutput', false);

boxplot2({prm(:, 1:maxNonEmptyBin)},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', ...
    'Actin Speed (nm/min)');

hAxes = get(hFig, 'CurrentAxes');
yTicks = 0:100:((ceil(max(prm(3,:) + prm(4,:)) / 100) + 1) * 100);
set(hAxes,'YTickLabel', num2str(yTicks'));
set(hAxes,'YTick',yTicks');
set(hAxes,'YLim', [0 yTicks(end) + 1]);

% Saving
fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_ActinSpeedVSdistanceToEdge.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%% Generate Adhesion length in fuction of distance to edge
data = num2cell(accLengthInDistanceBin, 1);
data = cellfun(@(c) nonzeros(c(~isnan(c))), data, 'UniformOutput', false);

prm = NaN(5, nDistanceBins);

bins = 1:nDistanceBins;
isNotEmptyBin = cellfun(@(c) numel(c) >= minSegmentsPerBin, data);
maxNonEmptyBin = max(bins(isNotEmptyBin));

for iBin = bins(isNotEmptyBin)
    sortedData = sort(data{iBin});
     
    prm(1,iBin) = sortedData(floor(numel(sortedData)/2)+1);
    prm(2,iBin) = sortedData(floor(numel(sortedData)/4)+1);
    prm(3,iBin) = sortedData(floor(3 * numel(sortedData)/4)+1);
    prm(4,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
    prm(5,iBin) = 1.5 * (prm(3,iBin) - prm(2,iBin));
     
    data{iBin} = sortedData;
end
 
hFig = figure('Visible', 'off');
hold on;

xlabels = arrayfun(@(iBin) [num2str(rangeDistance(iBin)) '-' num2str(rangeDistance(iBin+1))], ...
    1:maxNonEmptyBin, 'UniformOutput', false);

boxplot2({prm(:, 1:maxNonEmptyBin)},'color', [0.36 .63 .9], 'xlabels', xlabels, 'ylabel', ...
    'Adhesion Length (nm)');

hAxes = get(hFig, 'CurrentAxes');
yTicks = 0:500:((ceil(max(prm(3,:) + prm(4,:)) / 500)) * 500);
set(hAxes,'YTickLabel', num2str(yTicks'));
set(hAxes,'YTick',yTicks');
set(hAxes,'YLim', [0 yTicks(end) + 1]);

% Saving
fileName = fullfile(movieData.figures.directory, ...
    [getDirFromPath(movieData.imageDirectory) '_AdhLengthVSdistanceToEdge.eps']);
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%% END
movieData.figures.status = 1;
