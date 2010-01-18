function dataStruct = makiInitCoord(dataStruct, verbose, movieType, cutoffFix,timeRange)
%MAKIINITCOORD finds initial guesses for mammalian kinetochores
%
% SYNOPSIS: dataStruct = makiInitCoord(dataStruct)
%
% INPUT dataStruct: dataStruct as in makimakeDataStruct with at least the
%                       fields
%                         .rawMovieName - may contain the actual movie
%                         .rawMoviePath
%                         .dataProperties
%                   verbose (opt): 0: plot nothing
%                                  1: show progressText (default)
%                                  2: + show cutoff-plots
%                                  3: + show initial coordinate guesses
%                                     ordered by amplitude for every frame
%                   movieType: 1 for DV movies, 2 for STK movies.
%                              Optional. Default: 1.
%                   cutoffFix: 2-element vector indicating
%                           cutoff-type/cutoff-value for the selection of
%                           the "good" local maxima. Default: []
%                           This input is to be used for testing purposes
%                           only; a pre-set cutoff value should be stored
%                           in dataStruct.dataProperties.
%                           Setting [0,0] uses a pooled cutoff.
%                           KJ: I have modified makiMakeJob to allow user
%                           to enter cutoff, and it is now stored in
%                           dataStruct.dataProperties.detectionParam.
%                   timeRange: subset of timepoints for which makiInitCoord
%                           is being run. Optional. Default: [] (=all)
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
%                           This tends to be better than .amp
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
% created by: Jonas Dorn
% DATE: 28-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(verbose)
    verbose = 1;
end

if nargin < 3 || isempty(movieType)
    movieType = 1;
end

if nargin < 4 || isempty(cutoffFix)
    if isfield(dataStruct.dataProperties,'detectionParam')
        detectionParam = dataStruct.dataProperties.detectionParam;
        if isfield(detectionParam,'alphaValue')
            alphaValue = detectionParam.alphaValue;
            cutoffFix = [2 norminv(1-alphaValue)];
        else
            cutoffFix = [];
        end
    else
        cutoffFix = [];
    end
end

if nargin < 5 || isempty(timeRange)
    timeRange = [];
end

%=========================
%% COLLECT INPUT & SETUP
%=========================


if ischar(dataStruct.rawMovieName)
    rawMovie = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);
else
    rawMovie = [];
end
dataProperties = dataStruct.dataProperties;

% Jonas, 10/09: The background mask is only 1 sigma wide. This means that
% instead of Gaussian filtering, we were doing some kind of strange average
% filtering.

% setup filter parameters. 
backgroundFilterParms = dataProperties.FT_SIGMA * 15; % use *14 if you then filter the signalFiltered image to backgroundFilter to save some time
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

% get halfPatchSize to adjust centroid result. The center pixel of a 5x5x5
% array is (3,3,3), thus, we have to subtract that from centroid coordinates
halfPatchSize = dataProperties.FILTERPRM(4:6)/2+0.5;

% setup lists
nTimepoints = dataProperties.movieSize(4);
initCoordRaw = cell(nTimepoints,1); % raw initial coords (before testing)
% set up initCoord-Structure
tmp(1:nTimepoints,1) = struct('allCoord',[],...
    'allCoordPix',[],'correctionMu',[],'nSpots',[],...
    'initAmp',[],'amp',[],'data4MMF',[]);
dataStruct.initCoord = tmp;

goodTimes = 1:nTimepoints;
% check timeRange
if ~isempty(timeRange)
    goodTimes = intersect(goodTimes,timeRange);
end

minSpotsPerFrame = 20;


%=====================
%% MAIN LOOP
%=====================
if verbose
    progressText;
end
nTimes = numel(goodTimes); % not the best way to count, but it works.

%extract movie size from dataProperties
movieSize = dataProperties.movieSize;

for t=goodTimes(:)'
    
    % -- load images --
    % raw = current raw movie frame
    % filtered = current PSF-filtered movie frame
    % then calculate
    % background = current 10xPSF-filtered movie frame
    % amplitude = filtered - background
    % noise = locAvg(var(raw-filtered))
    
    if isempty(rawMovie)
        % movie has been passed directly
        raw = dataStruct.rawMovie(:,:,:,:,t);
    else
        switch movieType
            case 1 %DV files
                raw = cdLoadMovie({rawMovie,'raw'},[],struct('frames2load',{{t}},...
                    'crop',dataProperties.crop,'maxSize',dataProperties.maxSize));
            case {2,3} %STK files
                if movieType == 2
                    movieTPName = [rawMovie(1:end-5) num2str(t) '.STK']; %name of files storing current time point
                    imageStack = metaTiffRead(movieTPName,[],[],0); %read stack using Nedelec's code
                    raw = double(cat(3,imageStack.data)); %extract image information
                else %TIF series
                    movieTPName = [rawMovie(1:end-5) num2str(t) '.tif']; %name of files storing current time point
                    raw = []; %remove old image
                    for iSlice = 1 : movieSize(3) %read z-slices
                        raw(:,:,iSlice) = imread(movieTPName,iSlice);
                    end
                end
                cropInfo = dataProperties.crop; %cropping information
                if ~isempty(cropInfo)
                    cropInfo = cropInfo(:,1:2);
                    cropInfo(1,:) = max(cropInfo(1,:),[1 1]); %defaults for min
                    if cropInfo(2,1) == 0 %defaults for max
                        cropInfo(2,1) = size(raw,1);
                    end
                    if cropInfo(2,2) == 0
                        cropInfo(2,2) = size(raw,2);
                    end
                    raw = raw(cropInfo(1,1):cropInfo(2,1),cropInfo(1,2):cropInfo(2,2),:); %crop image
                end
        end
    end
    %remove offset
    offset = min(raw(:));
    raw = raw -  offset;
    % filter movie with gauss-filter
    
    % Jonas: use old version of border-correction to avoid risks (yes, I
    % have become paranoid about these things)
    filtered = fastGauss3D(raw,[],dataProperties.FILTERPRM(4:6),2,signalFilter);
    % filtering with sigma=15 is equal to filtering with 1 and then
    % with 14 if both filters are Gaussians. Thus, it should be possible to
    % make backgroundFilter smaller and filter the filtered instead of the
    % raw image.
    %background = fastGauss3D(filtered,[],backgroundFilterParms(4:6),2,backgroundFilter);
    background = fastGauss3D(raw,[],backgroundFilterParms(4:6),2,backgroundFilter);
    
    
    % amplitude is filtered image - background. This underestimates the
    % true signal. The underestimation becomes stronger if betterBackground
    % is not used.
    amplitude = filtered - background;
    
    % noise is local average (averaged over filter support) of squared
    % residuals of raw-filtered image
    noise = (raw-filtered).^2;
    noise = fastGauss3D(noise,[],dataProperties.FILTERPRM(4:6),2,noiseMask);
    
    % find local maxima
    locMax = loc_max3Df(filtered,[3 3 3]);
    locMaxIdx = sub2ind(dataProperties.movieSize(1:3),...
        locMax(:,1),locMax(:,2),locMax(:,3));
    
    % read amplitude and noise at local maxima
    initCoordTmp = [locMax, amplitude(locMaxIdx), noise(locMaxIdx), ...
        zeros(length(locMaxIdx),1)];
    initCoordTmp = sortrows(initCoordTmp,-4);
    % only take MAXSPOTS highest amplitudes. This will leave plenty of noise
    % spots, but it will avoid huge arrays
    initCoordTmp = initCoordTmp(1:min(dataProperties.MAXSPOTS,size(initCoordTmp,1)),:);
    
    %remove any spots with negative amplitude, because those are false
    %positives for sure - KJ
    if any(movieType == [2 3])
        initCoordTmp = initCoordTmp(initCoordTmp(:,4)>0,:);
    end
    
    % loop through all to get sub-pixel positions.
    raw = raw - background; %overwrite raw to save memory
    for iSpot = 1:size(initCoordTmp,1)
        
        % read volume around coordinate
        patch = stamp3d(...
            raw,dataProperties.FILTERPRM(4:6),...
            initCoordTmp(iSpot,1:3),0);
        
        %HLE,KJ - calculate low-index edge patch adjustment if relevant
        edgeAdjustTmp = initCoordTmp(iSpot,1:3) - halfPatchSize;
        edgeAdjustTmp = abs(edgeAdjustTmp) .* edgeAdjustTmp<0;
        
        % subpixel coord is integer coord plus centroid (subtract
        % halfPatchSize so that center coordinate of patch is (0,0,0))
        initCoordTmp(iSpot,1:3) = ...
            initCoordTmp(iSpot,1:3) + ...
            centroid3D(patch) - halfPatchSize;
        
        %HLE,KJ - correct position due to patch truncation due to proximity
        %to a low-index edge
        initCoordTmp(iSpot,1:3) = initCoordTmp(iSpot,1:3) + edgeAdjustTmp;
        
        % %         %KJ - Hunter's original attempt - has a bug:
        % %         %last line must have floor(edgeAdjustTmp)
        % %         %HLE - Adjust for edge effect. When patch is truncated by low-index
        % %         %edge of image, the halfPatchSize needs to be adjusted
        % %         edgeAdjustTmp = initCoordTmp(iSpot,1:3) - ...
        % %             (dataProperties.FILTERPRM(4:6)-1)./2 ;
        % %         initCoordTmp(iSpot,1:3) = initCoordTmp(iSpot,1:3) - ...
        % %             edgeAdjustTmp .* (edgeAdjustTmp < 0);
        
        % amplitude guess is integral.
        initCoordTmp(iSpot,6) = nanmean(patch(:));
        
    end
    
    % use maxPix-amplitude to calculate cutoff - meanInt is not consistent
    % with the rest of the measures!
    initCoordTmp = [initCoordTmp,...
        initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)),...
        initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)./max(initCoordTmp(:,4),eps))]; %#ok<AGROW>
    
    initCoordRaw{t} = initCoordTmp;
    
    if verbose
        progressText(find(t==goodTimes)/nTimes); % !! parfor counts backwards
    end
    
    if verbose > 2
        % plot frame with initial guesses for coords
        figure('Name',sprintf('frame %i',t))
        imshow(max(amplitude,[],3),[])
        hold on,
        for i=1:size(initCoordTmp,1),
            text(initCoordTmp(i,2),initCoordTmp(i,1),num2str(i),'Color','r');
        end
    end
end % loop timepoints



clear initCoordTmp
allCoord = cat(1,initCoordRaw{:});

%% Find cutoff

allCols = [4, 7, 8]; % columns used to calculate cutoff
% find cutoff based on amplitude/sqrt(noise/amp), though the others are
% very similar. Allow fallback if less than 25 spots per frame are found
% (this indicates that cutFirstHistmode of splitModes failed)
cutoff = zeros(3,1);
cutoff1 = splitModes(allCoord(:,4)); % amplitude
cutoff2 = splitModes(allCoord(:,7)); % amplitude/sqrt(nse) - dark noise
cutoff3 = splitModes(allCoord(:,8)); % amplitude/sqrt(nse/amp) - poisson
if ~isempty(cutoff1)
    cutoff(1) = cutoff1;
else
    cutoff(1) = NaN;
end
if ~isempty(cutoff2)
    cutoff(2) = cutoff2;
else
    cutoff(2) = NaN;
end
if ~isempty(cutoff3)
    cutoff(3) = cutoff3;
else
    cutoff(3) = NaN;
end


% plot all before selecting cutoff so that we can see what went wrong
if verbose > 1
    figure('Name',sprintf('cutoffs for %s',dataStruct.projectName))
    ah(1) = subplot(3,1,1);
    set(ah(1),'NextPlot','add')
    plot(ah(1),[1,nTimepoints],[cutoff(1) cutoff(1)]);
    xlabel('timepoints')
    ylabel('amplitude - signal in grey values')
    ah(2) = subplot(3,1,2);
    set(ah(2),'NextPlot','add')
    plot(ah(2),[1,nTimepoints],[cutoff(2) cutoff(2)]);
    xlabel('timepoints')
    ylabel('amplitude/sqrt(nse) - SNR for dark noise')
    ah(3) = subplot(3,1,3);
    set(ah(3),'NextPlot','add')
    plot(ah(3),[1,nTimepoints],[cutoff(3) cutoff(3)]);
    xlabel('timepoints')
    ylabel('amplitude/sqrt(nse/amp) - SNR for poisson noise')
    for t=goodTimes
        plot(ah(1),t,initCoordRaw{t}(:,4),'+')
        plot(ah(2),t,initCoordRaw{t}(:,7),'+')
        plot(ah(3),t,initCoordRaw{t}(:,8),'+')
    end
end

nn = nan(3,1);
% check for predetermined cutoff
if ~isempty(cutoffFix)
    % store old value
    oldCutoff = cutoff;
    % set cutoffIdx, cutoffCol
    cutoffIdx = cutoffFix(1);
    if cutoffIdx == 0
        cutoffCol = 0;
    else
        cutoffCol = cutoffIdx + 5;
        % update cutoff
        cutoff(cutoffIdx) = cutoffFix(2);
    end
else
    
    %     switch movieType
    %         case 1
    % note: ask only for spots in good frames
    
    if isfield(dataProperties,'minSpotsPerFrame')
        minN = dataProperties.minSpotsPerFrame;
    else
        minN = 20;
    end
    minGood = minN*length(goodTimes);
    nn(3) = sum(allCoord(:,8)>cutoff(3))/minGood;
    nn(2) = sum(allCoord(:,7)>cutoff(2))/minGood;
    nn(1) = sum(allCoord(:,4)>cutoff(1))/minGood;
    if nn(3) > 1
        cutoffIdx = 3;
        cutoffCol = 8;
    elseif nn(2) > 1
        cutoffIdx = 2;
        cutoffCol = 7;
    elseif nn(1) > 1
        cutoffIdx = 1;
        cutoffCol = 4;
    else
        error('less than %i spots per frame found. makiInitCoord failed',minSpotsPerFrame)
    end
    
    %         case {2,3}
    %
    %             %for the HMS data, the signal is very dim and photobleaches a lot,
    %             %and determining the cutoff on the fly does not work.
    %             %Thus, I'm fixing the cutoff criterion to the amplitude-to-noise
    %             %ratio and I'm fixing the cutoff to 2.32 (approximating a 0.01
    %             %alpha-value in a hypothesis test) - KJ
    %             %             cutoff(2) = 2.32;
    %             cutoff(2) = 1.645;
    %             %             cutoff(2) = 1.3;
    %             cutoffIdx = 2;
    %             cutoffCol = 7;
    %
    %     end
    
end

% remember the cutoff criterion used
dataStruct.statusHelp{3,3} = [cutoffIdx,cutoffCol];


if verbose > 1
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
% -- this turns out to be unnecessary due to alignment of coords



for t=goodTimes(:)'
    
    if cutoffCol > 0
        goodIdxL = initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx);
    else
        % do pooled cutoff. Accept if above at least one cutoff
        goodIdxL = bsxfun(@gt,initCoordRaw{t}(:,allCols),cutoff');
        goodIdxL = sum(goodIdxL,2)>1;
    end
    
    
    % count good spots
    dataStruct.initCoord(t).nSpots = sum(goodIdxL);
    
    
    
    % store pixel coords. Uncertainty is 0.25 pix
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
    if cutoffCol > 0
        twoAbove = find(initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx),2,'last');
        twoBelow = find(initCoordRaw{t}(:,cutoffCol) < cutoff(cutoffIdx),2,'first');
        dataStruct.initCoord(t).data4MMF = ...
            initCoordRaw{t}([twoAbove;twoBelow],1:3);
    end
end

% % code to determine optimal amplitude cutoff
% % check if we should just take average testStatistic instead of doing
% % complicated histogram stuff
% cordStruct = makiCoord2Cord(dataStruct.initCoord4MMF);
% [dataProperties] = ...
%     detectSpots_MMF_findAmplitudeCutoff(...
%     rawMovie, cordStruct, dataProperties,[],0)
