
function dataStruct = makiInitCoord(dataStruct, verbose)
%MAKIINITCOORD finds initial guesses for mammalian kinetochores
%
% SYNOPSIS: dataStruct = makiInitCoord(dataStruct)
%
% INPUT dataStruct: dataStruct as in makimakeDataStruct with at least the
%                   fields
%                     .rawMovieName
%                     .rawMoviePath
%                     .dataProperties
%                   verbose: whether or not to plot cutoff plots 
%                       (default: 0)
%
% OUTPUT dataStruct.initCoords: Structure of length nTimepoints with fields
%                       .allCoord     [x y z sx sy sz] coords and sigmas in
%                           microns, corrected for refocussing. Sigma is
%                           0.25 pixels
%                       .allCoordPix  same as coords, but in pixels and
%                           without correction.
%                       .correctionMu z-correction in microns
%                       .nSpots       number of spots
%                       .initAmp      maximum pixel intensity and estimated
%                           local noise of local maxima
%                       .amp          amplitude estimate from local
%                           integral around locMax
%                       .data4MMF     [x y z] pixel coordinates of two
%                           spots right above and right below the cutoff,
%                           which can be used with a function like
%                           detectSpots_MMF_findAmplitudeCutoff.m to
%                           calculate the t-test cutoff for mixture model
%                           fitting (see code snipped at end of fcn)
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

    % for conversion to microns in the end: get pixelSize
    pixelSize = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];

% get halfPatchSize to adjust centroid result. The center pixel of a 7x5
% array is (4,3), thus, we have to subtract that from centroid coordinates
halfPatchSize = dataProperties.FILTERPRM(4:6)/2+0.5;

% setup lists
nTimepoints = dataProperties.movieSize(4);
initCoordRaw = cell(nTimepoints,1); % raw initial coords (before testing)
% set up initCoord-Structure
tmp(1:nTimepoints,1) = struct('allCoord',[],...
    'allCoordPix',[],'correctionMu',[],'nSpots',[],...
    'initAmp',[],'amp',[],'data4MMF',[]);
dataStruct.initCoord = tmp;
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
    initCoordTmp = [locMax, amplitude(locMaxIdx), noise(locMaxIdx), ...
        zeros(length(locMaxIdx),1)];
    initCoordTmp = sortrows(initCoordTmp,-4);
    % only take 500 highest amplitudes. This will leave plenty of noise
    % spots, but it will avoid huge arrays
    initCoordTmp = initCoordTmp(1:min(dataProperties.MAXSPOTS,size(initCoordTmp,1)),:);
    
    
    % loop through all to get sub-pixel positions. 
    raw = raw - background; %overwrite raw to save memory
    for iSpot = 1:size(initCoordTmp,1)
        % read volume around coordinate
        patch = stamp3d(...
            raw,dataProperties.FILTERPRM(4:6),...
            initCoordTmp(iSpot,1:3),0);
        % subpixel coord is integer coord plus centroid (subtract
        % halfPatchSize so that center coordinate of patch is (0,0,0))
        initCoordTmp(iSpot,1:3)=...
            initCoordTmp(iSpot,1:3) + ...
            centroid3D(patch) - halfPatchSize;
        % amplitude guess is integral. 
        initCoordTmp(iSpot,6) = mean(patch(:));
    end
    
    % use maxPix-amplitude to calculate cutoff - meanInt is not consistent
    % with the rest of the measures!
     initCoordTmp = [initCoordTmp,...
            initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)),...
            initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)./max(initCoordTmp(:,4),eps))];
        initCoordRaw{t} = initCoordTmp;

    
    progressText(t/nTimepoints);
end % loop timepoints

clear initCoordTmp
allCoord = cat(1,initCoordRaw{:});


% find cutoff based on amplitude/sqrt(noise/amp), though the others are
% very similar. Allow fallback if less than 25 spots per frame are found
% (this indicates that cutFirstHistmode of splitModes failed)
cutoff = zeros(3,1);
cutoff(1) = splitModes(allCoord(:,4)); % amplitude
cutoff(2) = splitModes(allCoord(:,7)); % amplitude/sqrt(nse) - dark noise
cutoff(3) = splitModes(allCoord(:,8)); % amplitude/sqrt(nse/amp) - poisson
minGood = 25*nTimepoints;

if sum(allCoord(:,8)>cutoff(3)) > minGood
    cutoffIdx = 3;
    cutoffCol = 8;
elseif sum(allCoord(:,7)>cutoff(2)) > minGood
    cutoffIdx = 2;
    cutoffCol = 7;
elseif sum(allCoord(:,4)>cutoff(1)) > minGood
    cutoffIdx = 1;
    cutoffCol = 4;
else
    error('less than 25 spots per frame found. makiInitCoord failed')
end
% remember the cutoff criterion used
dataStruct.statusHelp{3,3} = [cutoffIdx,cutoffCol];

% plot all
if verbose == 2
    figure('Name',sprintf('cutoffs for %s',dataStruct.projectName))
    ah(1) = subplot(3,1,1);
    set(ah(1),'NextPlot','add')
    plot(ah(1),[1,nTimepoints],[cutoff(1) cutoff(1)]);
    ah(2) = subplot(3,1,2);
    set(ah(2),'NextPlot','add')
    plot(ah(2),[1,nTimepoints],[cutoff(2) cutoff(2)]);
    ah(3) = subplot(3,1,3);
    set(ah(3),'NextPlot','add')
    plot(ah(3),[1,nTimepoints],[cutoff(3) cutoff(3)]);
    for t=1:nTimepoints
        plot(ah(1),t,initCoordRaw{t}(:,4),'+')
        plot(ah(2),t,initCoordRaw{t}(:,7),'+')
        plot(ah(3),t,initCoordRaw{t}(:,8),'+')
    end
    title(ah(cutoffIdx),'cutoff selected')
    % if we want to look at this, we should do scatterCloud!!
%     figure('Name',sprintf('cutoff-comparison for %s',dataStruct.projectName))
%     minC = min(allCoord,[],1);
%     maxC = max(allCoord,[],1);
%     subplot(2,2,1)
%     plot(allCoord(:,4),allCoord(:,6),'.')
%     hold on, plot([minC(4),maxC(4)],[cutoff2,cutoff2])
%     hold on, plot([cutoff1, cutoff1],[minC(6),maxC(6)])
%     subplot(2,2,2)
%     plot(allCoord(:,7),allCoord(:,6),'.')
%     hold on, plot([minC(7),maxC(7)],[cutoff2,cutoff2])
%     hold on, plot([cutoff3, cutoff3],[minC(6),maxC(6)])
%     subplot(2,2,3)
%     plot(allCoord(:,4),allCoord(:,7),'.')
%     hold on, plot([minC(4),maxC(4)],[cutoff3,cutoff3])
%     hold on, plot([cutoff1, cutoff1],[minC(7),maxC(7)])
end



% loop and store only good locMax. Before that, get z-correction
% get correction from .log file
    % load $MOVIENAME.log, parse for start coordinate and determine focus
    % adjustment. This should give a z0-value for every frame that we can
    % use to correct the um-coords for tracking
    


for t=1:nTimepoints
    goodIdxL = initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx);
    
    % count good spots
    dataStruct.initCoord(t).nSpots = sum(goodIdxL);
    
    % store pixel coords. Uncertainty is 0.5 pix
    dataStruct.initCoord(t).allCoordPix = ...
        [initCoordRaw{t}(goodIdxL,1:3),...
        0.25*ones(dataStruct.initCoord(t).nSpots,3)];
    
    % store estimated amplitude and noise 
    dataStruct.initCoord(t).initAmp = initCoordRaw{t}(goodIdxL,4:5);
    
    % store integral amplitude
    dataStruct.initCoord(t).amp = [initCoordRaw{t}(goodIdxL,6),...
        zeros(dataStruct.initCoord(t).nSpots,1)];
    
    % store correction
    dataStruct.initCoord(t).correctionMu = 0; % for now
    
    % store coords in microns and correct
    dataStruct.initCoord(t).allCoord = ...
        dataStruct.initCoord(t).allCoordPix.*...
        repmat(pixelSize,dataStruct.initCoord(t).nSpots,2);
    dataStruct.initCoord(t).allCoord(:,3) = ...
        dataStruct.initCoord(t).allCoord(:,3) +...
        dataStruct.initCoord(t).correctionMu; 
    
    
    
    
    % store 2 spots above and below cutoff in case we want to get amplitude
    % cutoff for detector later. Note: the spots may be not exactly the two
    % above or below, since we sorted the data according to amplitudes
    % above.
    twoAbove = find(initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx),2,'last');
    twoBelow = find(initCoordRaw{t}(:,cutoffCol) < cutoff(cutoffIdx),2,'first');
    dataStruct.initCoord(t).data4MMF = ...
        initCoordRaw{t}([twoAbove;twoBelow],1:3);
end

% % code to determine optimal amplitude cutoff
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