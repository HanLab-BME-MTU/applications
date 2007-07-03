function dataStruct = makiInitCoord(dataStruct, verbose)
%MAKIINITCOORD finds initial guesses for mammalian kinetochores
%
% SYNOPSIS: dataStruct = makiInitCoord(dataStruct)
%
% INPUT dataStruct: dataStruct as in makimakeDataStruct with at least the
%                   fields
%                     .rawMovieName
%                     .rawMoviePath
%                     .filteredMovieName
%                     .dataProperties
%
% OUTPUT dataStruct.initCoords: t-by-1 cell array of n-by-5 arrays with
%                   initial guesses for coordinates (xyz), amplitudes and
%                   local noise. Contains only coordinates above the
%                   significance threshold
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 28-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(verbose)
    verbose = 0;
end

%=========================
%% COLLECT INPUT & SETUP
%=========================

rawMovie = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);
dataProperties = dataStruct.dataProperties;

% setup filter parameters
backgroundFilterParms = dataProperties.FT_SIGMA * 15;
backgroundFilterParms(4:6) = roundOddOrEven(backgroundFilterParms,'odd','inf');
% create background mask already
backgroundFilter = GaussMask3D(backgroundFilterParms(1:3),...
    backgroundFilterParms(4:6),[],1,[],[],1);
signalFilter = GaussMask3D(dataProperties.FILTERPRM(1:3),...
    dataProperties.FILTERPRM(4:6),[],1,[],[],1);

% make separated noise mask
noiseMask = {...
    ones(dataProperties.FILTERPRM(4),1,1)./dataProperties.FILTERPRM(4),...
    ones(1,dataProperties.FILTERPRM(5),1)./dataProperties.FILTERPRM(5),...
    ones(1,1,dataProperties.FILTERPRM(6))./dataProperties.FILTERPRM(6),...
    };

% setup lists
nTimepoints = dataProperties.movieSize(4);
initCoordRaw = cell(nTimepoints,1); % raw initial coords (before testing)
dataStruct.initCoord = cell(nTimepoints,1);

%=====================
%% MAIN LOOP
%=====================

progressText;

for t=1:nTimepoints

    % -- load images --
    % raw = current raw movie frame
    % filtered = current PSF-filtered movie frame
    % then calculate
    % background = current 10xPSF-filtered movie frame
    % amplitude = filtered - background
    % noise = locAvg(var(raw-filtered))

    raw = cdLoadMovie({rawMovie,'raw'},[],struct('frames2load',{{t}},...
        'crop',dataProperties.crop,'maxSize',dataProperties.maxSize));
    raw = raw - min(raw(:)); %remove offset
    % filter movie with gauss-filter
    filtered = fastGauss3D(raw,[],dataProperties.FILTERPRM(4:6),1,signalFilter);

    background = fastGauss3D(raw,[],backgroundFilterParms(4:6),1,backgroundFilter);
    
    amplitude = filtered - background;

    % noise is local average (averaged over filter support) of squared
    % residuals of raw-filtered image
    noise = (raw - filtered).^2;
    noise = fastGauss3D(noise,[],dataProperties.FILTERPRM(4:6),1,noiseMask);

    % find local maxima
    locMax = loc_max3Df(filtered,[3 3 3]);
    locMaxIdx = sub2ind(dataProperties.movieSize(1:3),...
        locMax(:,1),locMax(:,2),locMax(:,3));

    % read amplitude and noise at local maxima
    initCoordTmp = [locMax, amplitude(locMaxIdx), noise(locMaxIdx)];
    initCoordTmp = sortrows(initCoordTmp,-4);
    % only take 500 highest amplitudes. This will leave plenty of noise
    % spots, but it will avoid huge arrays
    initCoordTmp = initCoordTmp(1:dataProperties.MAXSPOTS,:);
    initCoordRaw{t} = [initCoordTmp,...
        initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)),...
        initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)./max(initCoordTmp(:,4),eps))];

    progressText(t/nTimepoints);
end % loop timepoints

clear initCoordTmp
allCoord = cat(1,initCoordRaw{:});

% find cutoff based on amplitude/sqrt(noise/amp), though the others are
% very similar
cutoff1 = splitModes(allCoord(:,4)); % amplitude
cutoff2 = splitModes(allCoord(:,6)); % amplitude/sqrt(nse) - dark noise
cutoff3 = splitModes(allCoord(:,7)); % amplitude/sqrt(nse/amp) - poisson

% plot all
if verbose == 2
    figure('Name',sprintf('cutoffs for %s',dataStruct.projectName))
    ah1 = subplot(3,1,1);
    set(ah1,'NextPlot','add')
    plot(ah1,[1,nTimepoints],[cutoff1 cutoff1]);
    ah2 = subplot(3,1,2);
    set(ah2,'NextPlot','add')
    plot(ah2,[1,nTimepoints],[cutoff2 cutoff2]);
    ah3 = subplot(3,1,3);
    set(ah3,'NextPlot','add')
    plot(ah3,[1,nTimepoints],[cutoff3 cutoff3]);
    for t=1:nTimepoints
        plot(ah1,t,initCoordRaw{t}(:,4),'+')
        plot(ah2,t,initCoordRaw{t}(:,6),'+')
        plot(ah3,t,initCoordRaw{t}(:,7),'+')
    end
    figure('Name',sprintf('cutoff-comparison for %s',dataStruct.projectName))
    minC = min(allCoord,[],1);
    maxC = max(allCoord,[],1);
    subplot(2,2,1)
    plot(allCoord(:,4),allCoord(:,6),'.')
    hold on, plot([minC(4),maxC(4)],[cutoff2,cutoff2])
    hold on, plot([cutoff1, cutoff1],[minC(6),maxC(6)])
    subplot(2,2,2)
    plot(allCoord(:,7),allCoord(:,6),'.')
    hold on, plot([minC(7),maxC(7)],[cutoff2,cutoff2])
    hold on, plot([cutoff3, cutoff3],[minC(6),maxC(6)])
    subplot(2,2,3)
    plot(allCoord(:,4),allCoord(:,7),'.')
    hold on, plot([minC(4),maxC(4)],[cutoff3,cutoff3])
    hold on, plot([cutoff1, cutoff1],[minC(7),maxC(7)])
end

% loop and store only good locMax
for t=1:nTimepoints
    goodIdx = initCoordRaw{t}(:,7) > cutoff3;
    dataStruct.initCoord{t} = initCoordRaw{t}(goodIdx,1:5);
    % store 2 spots above and below cutoff in case we want to get amplitude
    % cutoff for detector later. Note: the spots may be not exactly the two
    % above or below, since we sorted the data according to amplitudes
    % above.
    twoAbove = find(initCoordRaw{t}(:,7) > cutoff3,2,'last');
    twoBelow = find(initCoordRaw{t}(:,7) < cutoff3,2,'first');
    dataStruct.initCoord4MMF{t} = ...
        initCoordRaw{t}([twoAbove;twoBelow],1:4);
end

% % new code for detector
% % check if we should just take average testStatistic instead of doing
% % complicated histogram stuff
% cordStruct = makiCoord2Cord(dataStruct.initCoord4MMF);
% [dataProperties] = ...
%     detectSpots_MMF_findAmplitudeCutoff(...
%     rawMovie, cordStruct, dataProperties,[],0)








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutValue = splitModes(data)
% cut histogram where there is no data

% first guess via cutFirstHistMode
[cutIdx, cutVal,sp] = cutFirstHistMode(data,0);

% now check the local minima in the vicinity of the cutoff
spder = fnder(sp);
zeroList = fnzeros(spder);
zeroList = zeroList(1,:);
% evaluate
zeroVals = fnval(sp,zeroList);

% look in zeroList. Find one value before cutVal, three after. Go into
% zeroVals and find lowest minimum
[dummy,closestIdx] = min(abs(zeroList - cutVal));

% check only the minimas that are close by; two to the right and
% one to the left (don't forget that between two minima there will
% always be a maximum!)
indexList = (closestIdx-2):(closestIdx + 4);
indexList(indexList < 1 | indexList > length(zeroVals)) = [];

% also, never allow going left of the highest maximum!
[dummy,maxValIdx] = max(zeroVals);
indexList(indexList < maxValIdx) = [];

% find lowest
[dummy, cutIdx] = min(zeroVals(indexList));
% and determine break value
cutValue = zeroList(indexList(cutIdx));