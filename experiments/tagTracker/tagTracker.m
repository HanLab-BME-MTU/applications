function [idlisttrack,debugData] = tagTracker(movie,idlist,dataProperties,verbose,dbOpt)
%TAGTRACKER refines the linked tag positions found by the detector
%
% SYNOPSIS [idlisttrack] = tagTracker(movie,idlist,dataProperties,nSources)
%
% INPUT    movie: Can also be a cell with {moviename, movietype}
%          idlist
%          dataProperties
%
% OUTPUT   idlisttrack
%          debugData: structure with debug-fields
%
% tagTracker tries to map intensity distributions from source frames with
% separated tags to target frames with potentially non-separated tags to
% improve localization and resolution.
% It employs a Lucas-Kanade formalism.
%
% c: 2/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==================
%% TEST INPUT
%==================

% currently, all three input arguments are needed, and the dataProperties
% have not yet been updated - therefore, just test for nargin

% return an error if nargin does not equal 3 or 4
error(nargchk(3, 5, nargin, 'struct'));

% assign constants and parameters
constants.pixel2micron = [dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];
constants.filterSigma = dataProperties.FT_SIGMA;
% verbose: -1: nothing at all.
%           0: print to screen tracking/success
%           1: print start/end,
%           2: show parameter-convergence on screen
if nargin < 4 || isempty(verbose)
    constants.verbose = 0;
else
    constants.verbose = verbose;
end

% strategy options
constants.trackAll = 1;

% interpolation/objective function options
% interpolation {scheme, correctOrNot}
constants.interpolation = {'*cubic',1};
constants.gaussPixOnly = true;
constants.gradientOption = 1; % 1: filter, 2: don't filter
constants.correctBleaching = false;
constants.linearVariances = false; %whether or not to require a linear dependence of varXY on varZ from the detector

% debug options
% set defaults first
constants.numSources = 5;
debug = 0;
debugData = [];
if nargin < 5 || isempty(dbOpt)
    % do nothing
else
    debug = 1;
    if isstruct(dbOpt)
        % set number of sources (-to be moved into dataProperties)
        if isfield(dbOpt,'nSources')
            constants.numSources = dbOpt.nSources;
        end
        % set gradient option -> move into dataProperties
        if isfield(dbOpt,'gradientOption')
            constants.gradientOption = dbOpt.gradientOption;
        end
        if isfield(dbOpt,'gradientOptionQ')
            constants.gradientOptionQ = dbOpt.gradientOptionQ;
        end
        % data for the analysis of F-statistics:
        % see below
        if isfield(dbOpt,'fStats')
            debugData.fStats = [];
        end
        % tbd
        if isfield(dbOpt,'objectiveFunction')
            debugData.objectiveFunction = 1;
        end
        % track results:
        % see below
        if isfield(dbOpt,'trackResults')
            debugData.trackResults = [];
        end
        % fitting matrices, Qout, mse
        if isfield(dbOpt,'fitStats')
            debugData.fitStats = [];
        end
    end
end

% fitting options (tomlab is too slow)
% 'lsqnonlin'
% 'dom'
constants.fittingFcn = 'dom';

% calculate the  masks for taking and smoothening the gradient of the image
% differences
constants = gradients(constants,dataProperties);

% multiplicator of sigma for search radius
dataProperties.trackerRadiusMultiplicator = 25;
constants.trackerRadiusMultiplicator = ...
    dataProperties.trackerRadiusMultiplicator;


%=================


%=====================
%% INITIALIZE TRACKING
%=====================

% here, we read everything we need from idlist - we'll reconstruct it later
% (-> this makes it possible to use other input structures later)
% Then, we generate the tracking strategy and prepare the fitting matrices.

nTimepoints = length(idlist);

% remove first all single-occurence tags, so that we're only left with good
% ones. LG_deleteTag needs goodTimes, so we do a quick loop first
for t = nTimepoints:-1:1
    if ~isempty(idlist(t).linklist)
        goodIdx(t,1) = 1;
    end
end
goodTimes = find(goodIdx);

% since tags are good or bad over the entire movie, just check the first
% good time
tag2deleteIdx = find(ismember(idlist(goodTimes(1)).linklist(:,5),[2,3]));
% remove bad tags
if ~isempty(tag2deleteIdx)
    idlist = LG_deleteTag(idlist, tag2deleteIdx, goodTimes);
end

% now we can read out number of tags and assign vars
nTags = length(idlist(1).stats.labelcolor);

% assign variables.
nSpots = zeros(nTimepoints,1);
inputQmatrixDiag = zeros(nTimepoints,nTags,3);
inputCoords = zeros(nTimepoints,nTags,3);
inputAmp = zeros(nTimepoints,nTags);
p2m = repmat(constants.pixel2micron,[nTags,1]);

coords4Source = zeros(nTimepoints, 3*nTags);
searchRadius4Source = zeros(nTimepoints, 3*nTags);

% refit tag intensities to make sure we have the correct number of tags.
% Only overwrite intFit in the idlist!
idlist(1).stats.recalc = {2}; % only estimate
idTmp = linker(idlist,dataProperties);
idlist(1).stats.intFit = idTmp(1).stats.intFit;


for t = goodTimes'

    % get number of SPOTS
    nSpots(t) = max(idlist(t).linklist(:,2));

    % for qMatrix and coords: switch x and y!

    % multiply q by chi^2, which is stored in linklist(:,12). Error if
    % linklist isn't long enough - we can't use the old data

    % careful with all the reshaping!
    qMatDiag = reshape(reshape(...
        full(diag(idlist(t).info.detectQ_Pix)) .*...
        repeatEntries(idlist(t).linklist(:,12),3),...
        [3,nTags])',[1,nTags,3]);
    qMatDiag = qMatDiag(1,:,[2,1,3]);
    inputQmatrixDiag(t,:,:) = qMatDiag;
    % DONT USE SEARCHRADIUS ANYMORE
    % searchRadius: sqrt of Q * factor or from idlist.trackInit
    %     searchRadius = sqrt(reshape(squeeze(qMatDiag)',[3,nTags])') *...
    %         constants.trackerRadiusMultiplicator;
    %     if ~isempty(idlist(t).trackInit)
    %         searchRadius(idlist(t).trackInit(:,1),:) = ...
    %             idlist(t).trackInit(:,2:end);
    %     end
    %     searchRadius4Source(t,:) = reshape(searchRadius',[1,3*nTags]);


    % for estimated tags: use trackInit


    % Fill in coords for all frames: We'll need them as starting
    % positions or for search radii.

    % coords: convert to pixels and switch x,y
    coords = idlist(t).linklist(:,[10,9,11])./p2m; %nTags x 3
    inputCoords(t,:,:) = ...
        reshape(coords,[1,nTags,3]); %[1, nTags, 3]
    coords4Source(t,:) = reshape(coords',[1,3*nTags]);
    %     inputAmp(t,:,:) = reshape(idlist(t).linklist(:,8),[1,nTags,1]);
    % use "theoretical" intensities
    inputAmp(t,:,1) = idlist(1).stats.intFit.tagFactor .* ...
        exp(idlist(1).stats.intFit.xFit(end) * t);

end % loop nTimepoints

% find strategy. We can keep using the maximum number of spots to find
% sources (instead of dataProperties.MAXSPOTS), because all bad spots have
% been removed
stratInStruct = struct('nTimepoints',nTimepoints,...
    'nSpots',nSpots,'qVector',sum(sum(inputQmatrixDiag,3),2));
trackPairs = ...
    trackingStrategy2(stratInStruct,constants.trackAll,constants.numSources);

% build matrices for fitting
% We will solve the equation A*x=B where B has the covariance matrix V, and
% x are the positions of the tags. To have an reference position, the
% detected spot positions for source frames are also taken into account.
% Every tag will have its own matrices A,B and V, since there could be a
% case where both CENs are fused to their respective SPBs and only one of
% the two fusions can be resolved.
% V will be constructed as a vector, since we do not consider the
% cross-correlations between the tags in the fit (make sure it really plays
% no role!)
% Currently, I assume that the number of trackPairs will not change. If
% they do, the matrices have to be initialized as larger matrix to prevent
% slowing down of the code. Remaining zero-rows/cols in A are no problem
% for the code.

sourceList = unique(trackPairs(:,1));
nSources = length(sourceList);
nTrackPairs = size(trackPairs,1);

% A is the same for all three dimensions (it is also the same for all tags,
% but might not in the future
fittingMatrices(1:nTags) = struct(...
    'A',zeros(nSources+nTrackPairs,nTimepoints,1),...
    'B',zeros(nSources+nTrackPairs,1,3),...
    'V',zeros(nSources+nTrackPairs,1,3));

% fill in source data - this will produce zero-rows and cols in the fitting
% matrices. We will later remove them.
idx4A = sub2ind([nSources+nTrackPairs,nTimepoints],1:nSources,sourceList');
    for iTag = 1:nTags
        fittingMatrices(iTag).A(idx4A) = 1;
        fittingMatrices(iTag).B(1:nSources,1,:) = inputCoords(sourceList,iTag,:);
        fittingMatrices(iTag).V(1:nSources,1,:) = inputQmatrixDiag(sourceList,iTag,:);
    end


% fill sourceInfo
% sourceInfo (like targetInfo) is a structure containing all the
% information about the source/target necessary for the fitting process
% It is a nSources x nTags structure with fields
%   .centerCoord    : coordinates of the center
%   .coordList      : list of coordinates (xyz) belonging to the tag
%                       coordList-centerCoord gives the relative
%                       coordinate grid. If the source block is out of the
%                       image, it will be cut off. In the target image, the
%                       block will be wrapped around to conserve the number
%                       of residuals and to punish out-of-frame pixels
%                       Coordinates form a block around the
%                       center
%   .gaussList      : intensity of Gauss at coordinates in coordList
%   .ratioList      : relative intensity contribution of this tag to each
%                       "pixel" in coordList
%   .intList        : list of intensities interpolated from the images
%   .goodIdx        : in source: indices into coordList that belong to
%                       gauss, i.e. coordinates that
%                       really make up the source image. In target: indices
%                       into coordList, that both belong to gauss and are
%                       not wrapped
%   .amp            : amplitude of tag
%   .deltaImage     : intensity difference between source and target
%                       (target only)
%   .gradient       : gradient of difference image. Only calculated after
%                       fitting. Suffers from border effects, but the
%                       relevant residuals are well within the image
%   .size           : size of source image (source only)
%   .deltaCoord     : difference in coords from source to target (target
%                       only)
%   .background     : background of image. Calculated as mean of image with
%                       masked intensity distributions
%   .backgroundVar  : variance of background. Used to correct for the noise
%                       reduction caused by interpolation
%   .interpCorr     : correction to the residuals to offset noise reduction
%                       caused by interpolation. Only works for linear
%                       interpolation. (calculation is hard for cubic and
%                       probably impossible for spline)
%   .isSource       : 1 if it's a source


sourceInfo(1:nTimepoints,1:nTags) = struct(...
    'centerCoord',...
    mat2cell(coords4Source,ones(1,nTimepoints),3*ones(1,nTags)),...
    'searchRadius',...
    mat2cell(searchRadius4Source,ones(1,nTimepoints),3*ones(1,nTags)),...
    'coordList',[],'ratioList',[],'intList',[],'goodIdx',[],...
    'gaussList',[],'size',[],'deltaImage',[],...
    'amp',mat2cell(inputAmp,ones(1,nTimepoints),ones(1,nTags)),...
    'background',[],'backgroundVariance',[],...
    'deltaCoord',[], 'interpCorr',[],'isSource',0);
[sourceInfo(sourceList,:).isSource]=deal(1);

% loop and fill in source info
movieNoise = zeros(nTimepoints,1);
% remember randomState to keep analysis reproducible
randomState = rand('state');
rand('state',2756);
residualList = zeros(nTimepoints,nTags);
for t = sourceList'

    if iscell(movie)
        if strcmp(movie{2},'sedat')
            movieFrame = sedatLoadRaw(t,movie{1},dataProperties);
        else
            movieFrame = cdLoadMovie(movie,'',t);
        end
    else
        movieFrame = movie(:,:,:,:,t);
    end
    %disp(sprintf('preparing source %i',t))

    % estimate noise of target frame for subsequent F-statistics
    if movieNoise(t) == 0
        movieNoise(t) = imNoiseEstim(movieFrame);
        sourceInfo(t,1).backgroundVariance = movieNoise(t)^2;
    end

    sourceInfo(t,:) = extractIntensities(sourceInfo(t,:),[],movieFrame,constants);

    %---------- get the residuals of "bad" fits
    %
    % shift detected tags by 0.5*psf randomly along x,y,and z
    randomShift = rand(nTags*2,3);
    randomShift = (randomShift > 0.33) + (randomShift > 0.66) - 1;
    randomShift(all(randomShift == 0,2),:) = [];
    randomShift = randomShift(1:nTags,:);
    % shift by 3*gaussWidth which is approximately rayleigh
    randomShift = ...
        randomShift * 1 .* repmat(dataProperties.FT_SIGMA,nTags,1);
    % make parameter-vector
    randomShift = randomShift';
    randomShift = randomShift(:);

    % get residuals. no gradient
    [dummy,tmpInfo] = extractIntensities(sourceInfo(t,:), sourceInfo(t,:),...
        movieFrame, constants, randomShift, 0);

    for iTag = 1:nTags
        % if the position is outside the image, goodIdx is empty.
        % In that case, the residual is NaN.
        % Use sigmaResidual^2 instead of mean(residual)
        tmp = nansum(tmpInfo(iTag).deltaInt(tmpInfo(iTag).goodIdx).^2)/...
            (length(goodIdx)-nTags*3);
        residualList(t,iTag) = tmp;

    end
end
clear tmp % clear, b/c we want to make tmp into a structure below!
rand('state',randomState);

% the residuals will change with bleaching. Therefore, estimate the
% threshold using the intensity fit from the linker. Don't forget that
% we're looking at intensity^2! (background should not play a role here)
expValues = exp((1:nTimepoints)' * idlist(1).stats.intFit.xFit(end) * 2);
for iTag=1:nTags
    goodIdx = ~isnan(residualList(sourceList,iTag));
    expFactor = linearLeastMedianSquares(expValues(sourceList(goodIdx)),...
        residualList(sourceList(goodIdx),iTag),[],max(residualList(:,iTag)));
    % fill in residualList
    residualList(:,iTag) = expFactor * expValues;
end


%================================


%================================
%% TRACKING LOOP
%================================

% loop through the source/target pairs in trackingStrategy. For every
% target, loop until we have enough sources

% fittingIdx indicates the line where we will continue to write the fitting
% matrices
fittingIdx = nSources;
fittingIdx0 = fittingIdx; % 0th fitting idx


% DEBUG
if debug && isfield(debugData,'trackResults')

    % track results:
    % structure array (nTimepoints x nTags x numSources)
    % .startEndDelta 2x3 array with start delta, end delta
    % .sigma0 1x2 vector with start/end ssq of residuals
    % .sourcePos 1x3 array of source position
    % .sourceVar 1x3 array of variance of source position
    % .deltaVar  1x3 array of variance of delta
    % .info   ct, isSource, sourceT, targetT, exitflag, success
    % (in caller: add fields like deltaStartEndDelta as comparison to
    % reality)

    tmp(1:nTimepoints,1:nTags,1:constants.numSources) = ...
        struct('startEndDelta',[],'sigma0',[],...
        'sourcePos',[],'sourceVar',zeros(1,3),...
        'deltaVar',zeros(1,3),...
        'info',[]);
    debugData.trackResults = tmp;
    clear tmp

    % fill source pos, var
    [debugData.trackResults(sourceList,:,1).sourcePos] = ...
        deal(sourceInfo(sourceList,:).centerCoord);
    for t=sourceList'
        for j=1:nTags
            debugData.trackResults(t,j,1).sourceVar(:) = ...
                inputQmatrixDiag(t,j,:);
        end
    end
end
if debug && isfield(debugData,'fStats')

    % data for the analysis of F-statistics:
    % structure array (nTimepoints x nTags x numSources)
    % .residuals  vector of residuals
    % .resImGoodIdx
    % .initResiduals  vector of residuals of initial guess
    % .initResGoodIdx
    % .source     source image
    % .sourceImSize size of source image
    % .statistics  fprob, sigma0, targetNoise, n1, n2
    % .initStatistics  fprob, sigma0, targetNoise, n1, n2
    % .info   ct, isSource sourceT, targetT, exitflag, success

    tmp(1:nTimepoints,1:nTags,1:constants.numSources) = ...
        struct('residuals',[],'resImGoodIdx',[],'initResiduals',[],...
        'initResImGoodIdx',[],'source',[],'sourceImSize',[],...
        'sigma0',[],'statistics',[],...
        'initStatistics',[],'info',[]);
    debugData.fStats = tmp;
    clear tmp

    % fill source images
    [debugData.fStats(sourceList,:,1).source] = ...
        deal(sourceInfo(sourceList,:).intList);
    [debugData.fStats(sourceList,:,1).sourceImSize] = ...
        deal(sourceInfo(sourceList,:).size);

end

% initialize array for sigmaZero to estimate whether the residuals of the
% fit are acceptable. AmpRatio is purely for debugging purposes
% nFits x [s, t, sigmaResInit, sigmaRes, sigmaResInitEst, ampRatio] x nTags
collectedResiduals = zeros(nTimepoints*constants.numSources, 6, nTags);



% initialize counter for tracking (needed for debugging and display)
ct = 0;

movieFrameSize = prod(dataProperties.movieSize(1:3));

% loop through trackPairs. Use the ~isempty approach in case we switch to
% tracking until we find a match or something similar
while ~isempty(trackPairs)

    % select the part in the strategy which tracks the same target
    currentTarget = trackPairs(1,2);
    sameTargetIdx = find(trackPairs(:,2) == currentTarget);
    currentStrategy = trackPairs(sameTargetIdx,:);
    % ... and remove it from the global strategy
    trackPairs(sameTargetIdx,:) = [];

    % preassign targetInfo
    targetInfo = sourceInfo(currentTarget,:);
    targetIsSource = targetInfo(1).isSource;

    % load the movie frame for this target
    if iscell(movie)
        if strcmp(movie{2},'sedat')
            movieFrame = sedatLoadRaw(t,movie{1},dataProperties);
        else
            movieFrame = cdLoadMovie(movie,'',currentTarget);
        end
    else
        movieFrame = movie(:,:,:,:,currentTarget);
    end

    % estimate noise of target frame for subsequent F-statistics
    if movieNoise(currentTarget) == 0
        movieNoise(currentTarget) = imNoiseEstim(movieFrame);
        targetInfo(1).backgroundVariance = movieNoise(currentTarget)^2;
    end

    % prepare loop.
    successfulFits = 0;

    if debug
        % targetCt is needed for fStats and trackResults
        targetCt = 0;
    end

    % loop through currentStrategy until it is empty or we have found
    % enough sources
    while ~isempty(currentStrategy) && ...
            (constants.numSources - targetIsSource - successfulFits > 0)

        % count the track pair
        ct = ct + 1;


        % get current source
        currentSource = currentStrategy(1,1);
        if constants.verbose > -1
            trackingMsg = sprintf('tracking %i->%i (#%i/%i)',...
                currentSource,currentTarget,ct,nTrackPairs);
            disp(trackingMsg)
        end

        if debug
            % update targetCt
            targetCt = targetCt + 1;
        end

        % read initial parameters, and search radius
        % - search radius is not used right now

        [initialParameters, radius] = ...
            deal(zeros(3 * nTags,1));
        % initalParameters: Delta between acutal and estimated position
        % lowerBound: Delta - radius
        % upperBound: Delta + radius
        for iTag=1:nTags
            initialParameters((iTag-1)*3+1:iTag*3) = ...
                targetInfo(iTag).centerCoord - ...
                sourceInfo(currentSource,iTag).centerCoord;
            radius((iTag-1)*3+1:iTag*3) = ...
                targetInfo(iTag).searchRadius;
        end

        Jacobian = [];

        % initial sigmaResidual2 for each tag - calculate here in the
        % future
        [sigmaResidual2init,fRatioInit,fProbInit] = deal(zeros(nTags,1));

        % fit. Switch according to fitting function.
        switch constants.fittingFcn
            case 'lsqnonlin'
                % lsqnonlin

                %no bounds.
                lowerBound = [];
                upperBound = [];

                % set options
                if constants.verbose > 1
                    options = optimset('Jacobian','on',...
                        'Display','on','TolX',1E-3,'TolFun',1E-20,...
                        'TolPCG',0.01);
                else
                    options = optimset('Jacobian','on',...
                        'Display','off','TolX',1E-3,'TolFun',1E-20,...
                        'TolPCG',0.01);
                end


                % debug still:
                constants.movieNoise = movieNoise;
                constants.currentTarget = currentTarget;

                % generate anonymous function of parameters calling
                % track_lsqnonlinFitFcn
                optimFcn = @(parameters) (track_lsqnonlinFitFcn(...
                    parameters,sourceInfo(currentSource,:),targetInfo,movieFrame,constants));
                [parameters,resnorm,residuals,xfl,output,lambda,Jacobian] = ...
                    lsqnonlin(optimFcn ,initialParameters,...
                    lowerBound,upperBound,options);
                % set Jacobian to [] to use the unfiltered gradient
                % Jacobian = [];
                iter = -1;

            case 'dom'
                % use simple optimization which calculates parms
                % individually for each tag

                % loop settings
                maxIter = 100; % maximum number of iterations
                maxDelta = 0.01; % if the tags move less than this, quit

                % prepare fitting loop
                parameters = initialParameters;
                oldParameters = zeros(size(initialParameters,1),maxIter);
                iter = 0;
                xfl = 0; % exitflag. 0=no convergence

                %%%%%%%%%%%%%%%%%%%%%
                %  mostly debug
                %
                % get initial targetInfo
                [sourceInfo(currentSource,:), targetInfo] = ...
                    extractIntensities(sourceInfo(currentSource,:), targetInfo,...
                    movieFrame, constants, parameters, constants.gradientOption);

                % find ssq of residuals to later compare to tracking
                % result
                for iTag = 1:nTags
                    residuals = ...
                        targetInfo(iTag).deltaInt(targetInfo(iTag).goodIdx);
                    dofInit(iTag,1) = length(residuals)-length(parameters);
                    % since we did a subtraction of two noisy signals, the
                    % noise should increase by sqrt(2)
                    % potential improvement: calculate robust second moment
                    sigmaResidual2init(iTag,1) = (sum(residuals.^2)/(dofInit(iTag,1)));
                    fRatioInit(iTag,1) = (sigmaResidual2init(iTag)/2)/movieNoise(currentTarget)^2;
                    % no calculation of fprob for speed reasons
                    %fProbInit(iTag,1)  = fcdf(fRatioInit(iTag,1),dofInit(iTag,1),movieFrameSize);

                    %                     %f-test
                    %                     fRatio = sigmaResidual2(iTag,iter+1)/movieNoise(currentTarget)^2;
                    %                     fProb  = fcdf(fRatio,dof,movieFrameSize/2);
                    %                     isSuccess(iTag) = fProb < 0.95;

                    % debug
                    if debug && isfield(debugData,'fStats')

                        % store init - stuff
                        debugData.fStats(currentTarget,iTag,targetCt+targetIsSource).initResiduals = residuals;
                        debugData.fStats(currentTarget,iTag,targetCt+targetIsSource).initResImGoodIdx = targetInfo(iTag).goodIdx;
                        debugData.fStats(currentTarget,iTag,targetCt+targetIsSource).initStatistics = ...
                            [-1, sigmaResidual2init(iTag), movieNoise(currentTarget).^2,...
                            dofInit(iTag),movieFrameSize];

                    end

                end

                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % Loop until either 100 iterations, or no tag moved by more than 0.05 pix
                while iter < maxIter && xfl == 0;

                    % update iteration
                    iter = iter + 1;

                    % remember old parameters
                    oldParameters(:,iter) = parameters;

                    % get latest targetInfo
                    [sourceInfo(currentSource,:), targetInfo] = ...
                        extractIntensities(sourceInfo(currentSource,:), targetInfo,...
                        movieFrame, constants, parameters, constants.gradientOption);

                    % find new parameters for every single tag based on gradient and
                    % residuals
                    for iTag = 1:nTags

                        parameters( (iTag-1)*3+1:iTag*3 ) = ...
                            targetInfo(iTag).gradient \ targetInfo(iTag).deltaInt(...
                            targetInfo(iTag).goodIdx) ...
                            + parameters( (iTag-1)*3+1:iTag*3 );

                    end


                    % calculate how much parameters changed - if the tags moved by more
                    % than 0.01 pixel, go through another iteration

                    % To prevent oscillations from killing convergence,
                    % store all previous parameter combinations. If we come
                    % back to where we were before, stop.

                    % subtract current parameters from old parameters
                    delta = repmat(parameters,1,iter) - oldParameters(:,1:iter);
                    % reshape so that we get the 3d vector of each tag
                    delta = reshape(delta, 3, nTags, iter);
                    % calculate 3d distance and make delta into a
                    % iter-by-nTags array
                    delta = permute(sqrt(sum(delta.^2,1)),[2,3,1]);
                    % check for all distances below threshold
                    convergeIdx = find(all(delta < maxDelta, 1));

                    % if there is a convergeIdx, write into exitFlag how
                    % far back we had to go to find the same set of
                    % parameters. Else, continue loop.
                    % Run the optimization at least three times
                    if ~isempty(convergeIdx) && iter > 2
                        % there could be several convergeIdx. Take the last
                        % one.
                        xfl = iter - convergeIdx(end) + 1;
                    end
                end % Dom optimization loop


        end % switch tracking algorithm


        %                 trackerObjectiveFunction(...
        %                     movieFrame,sourceInfo(currentSource,:),targetInfo,...
        %                     initialParameters,constants,parameters,movieNoise(currentTarget));
        %                 set(gcf,'Name',sprintf('Optim %i',ct));
        %end

        %         if any(ct == [185,195,225,505])
        %         trackerObjectiveFunctionBig(...
        %             movieFrame,sourceInfo(currentSource,:),targetInfo,...
        %             initialParameters,constants,parameters,movieNoise(currentTarget));
        %         end


        % check output
        % I don't have a good idea on a good way to get individual
        % residuals for the tags, so I just calculate targetInfo again.
        % Also calculate gradient of difference image
        % !!!! get the jacobian that is actually used during calculation of
        % the parameters (=filtered)
        if isfield(constants,'gradientOptionQ')
            % allow setting which gradient is used by testing function.
            [dummy,targetInfo] = extractIntensities(sourceInfo(currentSource,:),...
                targetInfo,movieFrame,constants,parameters,constants.gradientOptionQ);
        else
            [dummy,targetInfo] = extractIntensities(sourceInfo(currentSource,:),...
                targetInfo,movieFrame,constants,parameters,constants.gradientOption);
        end

        % loop. If any one of the tags is bad, we do not keep the frame
        % Future: also check whether the lsqnonlin-optimization came to a
        % good result.
        % xfl: current exitFlag
        sigmaResidual2 = zeros(nTags,1);
        if xfl > 0

            dof = zeros(nTags,1);
            isSuccess = 1;
            iTag = 0;
            %             while isSuccess && iTag < nTags
            for iTag = 1:nTags

                % use only residuals of Gauss!
                residuals = ...
                    targetInfo(iTag).deltaInt(targetInfo(iTag).goodIdx);
                dof(iTag) = length(residuals)-length(parameters);
                % since we did a subtraction of two noisy signals, the
                % noise should increase by sqrt(2)
                % potential improvement: calculate robust second moment
                sigmaResidual2(iTag) = (sum(residuals.^2)/(dof(iTag)));

                % f-test
                fRatio(iTag,1) = (sigmaResidual2(iTag)/2)/movieNoise(currentTarget)^2;
                fProb(iTag,1)  = fcdf(fRatio(iTag),dof(iTag),movieFrameSize); %removed /2 for n-noise
                isSuccess = (fProb(iTag) < 0.95) && isSuccess;

                % debug
                if debug && isfield(debugData,'fStats')

                    debugData.fStats(currentTarget,iTag,targetCt+targetIsSource).residuals = residuals;
                    debugData.fStats(currentTarget,iTag,targetCt+targetIsSource).resImGoodIdx = targetInfo(iTag).goodIdx;
                    debugData.fStats(currentTarget,iTag,targetCt+targetIsSource).statistics = ...
                        [fProb(iTag), sigmaResidual2(iTag), movieNoise(currentTarget).^2,...
                        dof(iTag),movieFrameSize];

                end


            end



        else
            isSuccess = 0;
            fProb = -1;
        end





        % !!!!!!!!!!!!!!!!!!!!

        % set success to 1 if the sum of sigmaResidual2 has been lowered
        if (sum(sigmaResidual2) < sum(sigmaResidual2init)) && all(fProb > 0)
            isSuccess = 1;
        else
            isSuccess = 0;
        end

        % debug - store .info
        if debug
            if isfield(debugData,'trackResults')
                debugData.trackResults(currentTarget,1,targetCt+targetIsSource).info = ...
                    [ct, targetIsSource, currentSource,currentTarget, xfl, iter, isSuccess];
            end
            if isfield(debugData,'fStats')
                debugData.fStats(currentTarget,1,targetCt+targetIsSource).info = ...
                    [ct, targetIsSource, currentSource,currentTarget, xfl, iter, isSuccess];
            end
        end


        if constants.verbose > 0
            disp(sprintf(['init =',repmat(' %1.4f',1,length(parameters(:))),...
                ' res =',repmat(' %1.4E',1,nTags),...
                ' prob =',repmat(' %1.4f',1,nTags)],...
                initialParameters'+cat(2,sourceInfo(currentSource,:).centerCoord),...
                sigmaResidual2init',...
                fProbInit'))
            disp(sprintf(['end  =',repmat(' %1.4f',1,length(parameters(:))),...
                ' res =',repmat(' %1.4E',1,nTags),...
                ' prob =',repmat(' %1.4f',1,nTags)],...
                parameters'+cat(2,sourceInfo(currentSource,:).centerCoord),...
                sigmaResidual2',...
                fProb'))
        end
        if constants.verbose > -1
            ans ={'no success', 'success'};
            disp(sprintf('%s, %i (%i), %1.4f\n',ans{isSuccess+1},xfl,iter, max(fProb)))
        end

        % if the tracking was a success, we go on
        if isSuccess

            % count the fit
            successfulFits = successfulFits + 1;

            % update the index here, since there is a fittingMat for each
            % tag
            fittingIdx = fittingIdx + 1;
            % loop through tags and write the fittingMatrices. Multiply V
            % by sigmaResidual
            for iTag = 1:nTags
                % parameters = target - source
                %fittingMatrices(iTag).A(fittingIdx,currentTarget,1) = 1;
                %fittingMatrices(iTag).A(fittingIdx,currentSource,1) = -1;
                fittingMatrices(iTag).A(fittingIdx,[currentSource,currentTarget],1) = [-1,1];
                fittingMatrices(iTag).B(fittingIdx,1,:) = ...
                    targetInfo(iTag).deltaCoord(:);
                % if we come from lsqnonlin: use its Jacobian
                if isempty(Jacobian)
                    qMat = inv(...
                        targetInfo(iTag).gradient'...
                        * targetInfo(iTag).gradient);
                else
                    qMat = inv(...
                        Jacobian(:,((iTag-1)*3+1):(iTag*3))'...
                        * Jacobian(:,((iTag-1)*3+1):(iTag*3)));
                end
                qMatDiag = diag(qMat) * sigmaResidual2(iTag);
                fittingMatrices(iTag).V(fittingIdx,1,:) = qMatDiag(:);
            end

            % store sigmaResidual2, sigmaResidual2Init for later estimation
            % of good fits
            % careful: collectIdx is used below to discard excess rows
            collectIdx = fittingIdx - fittingIdx0;
            collectedResiduals(collectIdx, 1, :) = currentSource;
            collectedResiduals(collectIdx, 2, :) = currentTarget;
            collectedResiduals(collectIdx, 3, :) = ...
                reshape(sigmaResidual2init,[1,1,nTags]);
            collectedResiduals(collectIdx, 4, :) = ...
                reshape(sigmaResidual2,[1,1,nTags]);
            % for debugging purposes, store ampRatio
            if debug
                for iTag = 1:nTags
                    collectedResiduals(collectIdx, 6, iTag) = ...
                        targetInfo(1,iTag).amp/...
                        sourceInfo(currentSource,iTag).amp;
                end
            end

        end %if isSuccess


        % debug - fill trackResults here, b/c of Q
        if debug && isfield(debugData,'trackResults')
            for iTag = 1:nTags
                debugData.trackResults(currentTarget,iTag,targetCt+targetIsSource).startEndDelta =...
                    [initialParameters((iTag-1)*3+1:iTag*3)';...
                    parameters((iTag-1)*3+1:iTag*3)'];
                debugData.trackResults(currentTarget,iTag,targetCt+targetIsSource).sigma0(1) = ...
                    sigmaResidual2init(iTag);
                debugData.trackResults(currentTarget,iTag,targetCt+targetIsSource).sigma0(2) = ...
                    sigmaResidual2(iTag);

                if isSuccess

                    debugData.trackResults(currentTarget,iTag,targetCt+targetIsSource).deltaVar(:) = ...
                        fittingMatrices(iTag).V(fittingIdx,1,:);
                end
            end

        end


        % whether it was successful or not, we remove the current s/t pair
        % from the strategy

        currentStrategy(1,:) = [];

    end % loop currentStrategy
    
%     % debug
%     if mod(collectIdx,10)==0
%         save(sprintf('fittingMatrices_A_%i_%s.mat',collectIdx,nowString),'fittingMatrices');
%     end


end % while ~isempty trackingStrategy

%====================================

%save(sprintf('fittingMatrices_A_%s.mat',nowString),'fittingMatrices');

%====================================
%% REFINE TRACKING RESULTS
%====================================

% First, check for good trackings: Are the residuals reasonably low?
% (Fit collected initial residuals, then remove all whose residuals are
% higher than expected)
% Secondly, transform the precision so that the measurements from the
% detector have weight 1 - theoretically, this means to divide Var(xy) by
% Mean(Var(xy)), and the same for z; instead we just set all to 1. Take the
% robust mean and kick out the observations from the detector that are
% considered to be outliers. Then, divide the Variances from the tracker by
% the mean precision from the detector. This will set sigma0(detector) to
% pixel, and it will give a meaning to the output MSE.

% 1) estimate the expected sum of residuals

% remove excess entries
collectedResiduals(collectIdx+1:end,:,:) = [];

% read sigmaResidual2Init
Y = sum(collectedResiduals(:,3,:),3);
% construct A. Put ones wherever sigma corresponds to that tag
%oneZero = (mat2cell(ones(collectIdx,nTags),collectIdx,ones(1,nTags)));
A = [ones(size(Y)),collectedResiduals(:,2,1)];
[xFit,sigmaX] = robustExponentialFit2(Y,A,constants.verbose > 0);


% estimatedResiduals are A*X
estimatedResiduals = exp(A * [log(xFit(1));xFit(end)]);
%collectedResiduals(:,5,1) = reshape(estimatedResiduals,[collectIdx,1,nTags]);

% if sigmaResidualEnd is above the estimatedResiduals + 10%, don't count
% (this will remove quite a number of s-s trackings. However, these are
% probably not very good, anyway)
highResidualIdx = find(sum(collectedResiduals(:,4,:),3) > estimatedResiduals * 1.1);

if constants.verbose > 0
    % plot and save
    set(gcf,'Name','total ssq of residuals vs. target frame');
    hold on, plot(A(:,end),estimatedResiduals * 1.1, 'r')
    hold on,plot(A(highResidualIdx,end),Y(highResidualIdx),'og'),
    fh = gcf;
    saveas(fh,sprintf('ssqVsTarget_%s.fig',nowString))
end
highResidualIdx = highResidualIdx + fittingIdx0;


% remove rows below (where we also remove other rows)


% 2) transform precision - only if more than 3 sources overall

if nSources > 3
% find mena detectorPrecision in XY and in Z. Detect very high
% uncertainties
% !!! fit precision with exponential - uncertainties increase with
% decreasing intensity. Remove outliers
detectorPrecisionXY = catStruct(4,'fittingMatrices.V(:,:,1:2)');
detectorPrecisionXY = squeeze(detectorPrecisionXY(1:fittingIdx0,:,:,:));
detectorPrecisionZ = catStruct(4,'fittingMatrices.V(:,:,3)');
detectorPrecisionZ = squeeze(detectorPrecisionZ(1:fittingIdx0,:,:,:));
% dp* is now a nTimepoints by nDims by nTags array. Since we're fitting
% for each tag individually, we will calculate an estimated precision for
% each tag.
for iTag=nTags:-1:1 % so that we don't need to preallocate
    % varXY and varZ show a very strong linear dependence. Therefore, we
    % can fit both with the same exponential.
    varX = detectorPrecisionXY(:,1,iTag);
    varY = detectorPrecisionXY(:,2,iTag);
    varZ = detectorPrecisionZ(:,iTag);
    t = (1:fittingIdx0)'; % t is equivalent to row for the detected spots

    % remove zero - entries
    zeroIdx = ~varX;
    varX(zeroIdx) = [];
    varY(zeroIdx) = [];
    varZ(zeroIdx) = [];
    t(zeroIdx) = [];

    % find frames where the dependence between varXY and varZ is linear
    % Problem: the linear dependence is so excellent that there are a lot
    % of frames removed

    if constants.linearVariances
        [x,dummy,goodRows] = linearLeastMedianSquares([[varX;varY],ones(2*length(t),1)],[varZ;varZ]);
        aa=[[varX;varY],ones(2*length(t),1)];bb=[varZ;varZ];
        badIdx = setdiff(1:2*length(t),goodRows);
        % remove bad data
        rejectIdx=unique((mod(badIdx-1,length(t))+1));

        if constants.verbose > 0
            fh = figure('Name',sprintf('Variances %s',idlist(1).stats.labelcolor{iTag}));
            ah=subplot(2,1,1);
            plot(aa(:,1),bb,'.')
            hold on, plot(aa(:,1),aa*x,'r')

            hold on, plot(aa(badIdx,1),bb(badIdx),'or')
            set(get(ah,'Title'),'String',sprintf('%f ',x))
            saveas(fh,sprintf('DetectorVariance_tag%i_%s.fig',iTag,nowString));
            subplot(2,1,2)
            plot(repmat(t,2,1),[varX;varY],'.b',t,varZ,'.g')
            % rejected: red circles
            hold on, plot(repmat(t(rejectIdx),3,1),...
                [varX(rejectIdx);varY(rejectIdx);varZ(rejectIdx)],'or')
            % cyan/magenta dots for source/target
            plot(sourceList,zeros(size(sourceList)),'.c')
        end

        varX(rejectIdx) = [];
        varY(rejectIdx) = [];
        varZ(rejectIdx) = [];
        t(rejectIdx) = [];


    else
        % don't do anything
    end

    % exponential fit. t is just a row-index; time is defined by sourceList
    A = [blkdiag(ones(length(t)*2,1),ones(size(t))),repmat(sourceList,3,1)];
    B = [varX;varY;varZ];
    % the exponential fit relies on fminsearch, which can fail, so run it
    % with only xy first, and then use the result as an initial guess for
    % the constrained fit. 
    x0 = robustExponentialFit2([varX;varY],A(1:length(t)*2,[1,3]),0);
    x0 = x0([1,1,2]);
    [xFit, dummy, goodRows] = robustExponentialFit2(B,A,constants.verbose>0,x0);
    if constants.verbose>0
        set(gcf,'Name',sprintf('detector variance xy,z of tag %i',iTag))
    end

    % get goodTime
    varianceEstimators{iTag} = xFit;
    % careful: t is not timepoints, but fit-indices of the source-rows
    % goodRows is e.g. [1:30,1:30,1:30]. Therefore, use mod, assuming that
    % if one dim is good, all are.
    goodDetection{iTag} = t(unique((mod(goodRows-1,length(t))+1)));
    badDetection{iTag} = setdiff(t,goodDetection{iTag});


    % write estimated errors. Divide by the estimated detection error at
    % t=1 for both detected and tracked tags. This sets weight 1 for the
    % detection in the first frame.
    fittingMatrices(iTag).V(goodDetection{iTag},1,1) = ...
        xFit(1) * exp(xFit(end) * goodDetection{iTag});
    fittingMatrices(iTag).V(goodDetection{iTag},1,2) = ...
        xFit(1) * exp(xFit(end) * goodDetection{iTag});
    fittingMatrices(iTag).V(goodDetection{iTag},1,3) = ...
        xFit(2) * exp(xFit(end) * goodDetection{iTag});

    fittingMatrices(iTag).V(:,1,1:2) = ...
        fittingMatrices(iTag).V(:,1,1:2) ./ xFit(1) * exp(xFit(end));
    fittingMatrices(iTag).V(:,1,3) = ...
        fittingMatrices(iTag).V(:,1,3) ./ xFit(2) * exp(xFit(end));

end

else
   % if only 3 or fewer sources, we still need to fill in bad 
   [badDetection{1:nTags}] = deal([]);
   varianceEstimators = [];
end

%====================================



%====================================
%% CALCULATE POSITIONS
%====================================

% at this point, we basically just have to run myLscov to get the
% xyz-coords of the tags in the target frames.
% before, we need to remove the empty cols in A, though. These will make up
% the list of frames that could not be tracked.



% because of the way the new lscov is coded, it cannot cope with zero-rows
% anymore -> remove them. If we're removing detected spots, then goodRows
% can be different for different tags
failIdx = (all(~fittingMatrices(1).A(:,:,1),2));
clear goodRows
[goodRows{1:nTags}] = deal(1:size(fittingMatrices(1).A,1));
for iTag = 1:nTags
    goodRows{iTag}([find(failIdx); highResidualIdx; badDetection{iTag}]) = [];
end

% since all A-matrices are the same, we can find the failIdx on a single
% one
clear goodIdx
[goodIdx{1:nTags}] = deal((1:nTimepoints)');
for iTag = 1:nTags
    failIdx = find(all(~fittingMatrices(1).A(goodRows{iTag},:,1),1));
    goodIdx{iTag}(failIdx) = [];
end

% prepare "output"
outputCoords = zeros(nTimepoints,nTags,3);
[outputQmatrixDiag,mse] = deal(zeros(nTimepoints,3,nTags));


% loop through dimensions and fit
for iTag = 1:nTags
    for dim=1:3
        [outputCoords(goodIdx{iTag},iTag,dim), outputQmatrixDiag(goodIdx{iTag},dim,iTag), mse(goodIdx{iTag},dim,iTag)] = ...
            myLscov(fittingMatrices(iTag).A(goodRows{iTag},goodIdx{iTag},1),...
            fittingMatrices(iTag).B(goodRows{iTag},1,dim),fittingMatrices(iTag).V(goodRows{iTag},1,dim));
    end
end

if constants.verbose > -1
    if isempty(varianceEstimators)
        % can't write the tilde on this keyboard...
    else
    fh=figure('Name',sprintf('MSE %s x;y;z',sprintf('%s ',idlist(1).stats.labelcolor{:})));
    for dim=1:3,
        for iTag=1:nTags,
            ah=subplot(3,nTags,(dim-1)*nTags+iTag);
            plot(mse(:,dim,iTag),'.'),
            set(get(ah,'Title'),'String',sprintf('s0 %f',varianceEstimators{iTag}(floor(dim/3)+1)));
        end,
    end
    saveas(fh,sprintf('mse_%s.fig',nowString))
    end
end

if constants.verbose > 0
    fh=figure('Name',sprintf('Qout %s x;y;z',sprintf('%s ',idlist(1).stats.labelcolor{:})));
    for dim=1:3,
        for iTag=1:nTags,
            subplot(3,nTags,(dim-1)*nTags+iTag);
            plot(goodIdx{iTag},outputQmatrixDiag(goodIdx{iTag},dim,iTag).^2,'.'),

        end,
    end
    saveas(fh,'varOut.fig');

    fh=figure('Name',sprintf('V input %s x;y;z',sprintf('%s ',idlist(1).stats.labelcolor{:})));
    for dim=1:3,
        for iTag=1:nTags,subplot(3,nTags,(dim-1)*nTags+iTag),
            plot(fittingMatrices(iTag).V(:,1,dim),'.')
            hold on, plot(goodRows{iTag},fittingMatrices(iTag).V(goodRows{iTag},1,dim),'or'),
        end,
    end
    saveas(fh,'varIn.fig');
end
% !!! Adjust the variance matrix. For whatever reason, the input error
% estimate is better than the output. In other words, we claim that the
% sigma zero is better estimated by detector and tracker than by the
% fitting, so we divide the output variance by sigma zero (mse).
for iTag = 1:nTags
    outputQmatrixDiag(goodIdx{iTag},:,iTag) = ...
        outputQmatrixDiag(goodIdx{iTag},:,iTag).^2; % removed division by mse 7/14
end
% figure('Name',sprintf('MSE %s x;y;z',sprintf('%s ',idlist(1).stats.labelcolor{:})));
% for dim=1:3,for iTag=1:nTags,subplot(3,nTags,(dim-1)*nTags+iTag),plot(outputQmatrixDiag(:,dim,iTag),'.'),end,end

% figure,plot(outputQmatrixDiag);
% hold on, plot(squeeze(inputQmatrixDiag),'.')
% sa=squeeze(sum(abs(fittingMatrices.A),1));
% hold on, plot(outputQmatrixDiag.*sqrt(sa),'o' );
% for iTag = 1:nTags
%     sa=repmat((sum(abs(fittingMatrices(iTag).A),1))',[1,3]);
%     outputQmatrixDiag(:,:,iTag) = outputQmatrixDiag(:,:,iTag) .* (sa);
% end

% DEBUG
%save(sprintf('trackerWorkspace%s.mat',nowString));

% if debug: shorten the matrices so that it isn't insanely big
if debug && isfield(debugData,'fitStats')
    debugData.fitStats.fittingMatrices = fittingMatrices;
    for iTag=1:nTags
        debugData.fitStats.fittingMatrices(iTag).A = ...
            (fittingMatrices(iTag).A(goodRows,:,1));
        debugData.fitStats.fittingMatrices(iTag).A = ...
            sparse(fittingMatrices(iTag).A);
        debugData.fitStats.fittingMatrices(iTag).B = ...
            fittingMatrices(iTag).B(goodRows,:,:);
        debugData.fitStats.fittingMatrices(iTag).V = ...
            fittingMatrices(iTag).V(goodRows,:,:);
    end
    % all A-matrices are identical
    for iTag=2:nTags
        debugData.fitStats.fittingMatrices(iTag).A = [];
    end

    debugData.fitStats.mse = mse;
    debugData.fitStats.outputQmatrixDiags = outputQmatrixDiag;
    debugData.fitStats.goodIdx = goodIdx;
end

%===================================

%===================================
%% WRITE OUT IDLISTTRACK
%===================================

% loop through idlist. Fill in the new coords and q-matrices. Keep the old
% entry wherever we had a lasting problem

idlisttrack = idlist;

% merge goodIdx into a single list of indices with frames where all
% tracking succeeded
goodIdxAll = (1:nTimepoints)';
for iTag=1:nTags
    goodIdxAll = intersect(goodIdxAll,goodIdx{iTag});
end

for t = [goodTimes';NaN,goodTimes(1:end-1)']
    if any(t(1)==goodIdxAll)

        % write new coords
        newCoords = outputCoords(t(1),:,:);

        % theoretically, the tags are still correctly sorted. However, it
        % is possible that there are swaps for tags with very similar
        % intensities. Should need arise, we could simply associate the new
        % tags with the old ones calling LAP (don't do that for the moment,
        % though).

        % newCoords is 1x2x3 - just squeeze and multiply
        % ... and don't forget to swap x and y again
        idlisttrack(t(1)).linklist(:,[10,9,11]) = ...
            permute(newCoords,[2,3,1]).*p2m;

        % write new Q-matrices. Divide the covariance matrix by the old
        % sigma0, so that we get the correct covariance once we multiply,
        % and so that we conserve the sigma0 to get the correct intensities
        newQdiag = squeeze(outputQmatrixDiag(t(1),:,:))./repmat(idlisttrack(t(1)).linklist(:,12),[1,3])';
        idlisttrack(t(1)).info.totalQ_Pix = diag(newQdiag(:));

        % write sigma0 of the tracker. Since it differs depending on the
        % dimension, we have to add cols 13:15
        idlisttrack(t(1)).linklist(:,13:15) = squeeze(mse(t(1),:,:))';

        % write list of sources
        targetIdx=fittingMatrices(1).A(:,t(1),1) == 1;
        [dummy, sourceList]=find(fittingMatrices(1).A(targetIdx,:)==-1);
        idlisttrack(t(1)).info.sourceList = sourceList;

        % write new spotNum
        idlisttrack(t(1)).linklist(:,2) = (1:nTags)';

        % update flags: Wherever we had a 1 in ll-3, there should be a 2
        % now, since we tracked all the estimated spots
        estimateIdx = idlisttrack(t(1)).linklist(:,3) == 1;
        idlisttrack(t(1)).linklist(estimateIdx,3) = 2;

    else

        % fill linklist 13-15 with zeros if necessary
        idlisttrack(t(1)).linklist(:,13:15) = zeros(nTags,3);

    end


    % write new linkup/linkdown. Do this for all frames!
    if ~isnan(t(2))
        idlisttrack(t(2)).linklist(:,7) = idlisttrack(t(1)).linklist(:,2);
        idlisttrack(t(1)).linklist(:,6) = idlisttrack(t(2)).linklist(:,2);
    end


end



%=========================================================================

%% Subfunction gradients
function constants = gradients(constants,dataProperties)


% GaussGrad-kernel. Use 1/2 sigma to filter - the difference between the
% images leads to much higher image frequencies than just the single psfs!
halfSigmaXZ = dataProperties.FILTERPRM([1,3])/2;
filterSizeXZ = roundOddOrEven(4*halfSigmaXZ,'odd','inf');

% calculate four kernels for faster calculation of the derivative later

x=(-filterSizeXZ(1)/2:filterSizeXZ(1)/2)./halfSigmaXZ(1);
z=(-filterSizeXZ(2)/2:filterSizeXZ(2)/2)./halfSigmaXZ(2);
% erf is integral of Gauss
constants.GaussXY = diff(0.5 * erfc(-x./sqrt(2)));
constants.GaussZ = diff(0.5 * erfc(-z./sqrt(2)));
% for gradient: Gauss is integral of GaussGrad
constants.GaussGradXY = diff(exp(-x.^2/2))/(sqrt(2*pi)*halfSigmaXZ(1));
constants.GaussGradZ = diff(exp(-z.^2/2))/(sqrt(2*pi)*halfSigmaXZ(2));

% reshape Z
constants.GaussZ = reshape(constants.GaussZ,1,1,filterSizeXZ(2));
constants.GaussGradZ = reshape(constants.GaussGradZ,1,1,filterSizeXZ(2));

constants.GaussSigma = halfSigmaXZ([1,1,2]);
constants.GaussSize = filterSizeXZ([1,1,2]);



