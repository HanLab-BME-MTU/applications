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
constants.gaussPixOnly = 1;
constants.gradientOption = 1; % 1: filter, 2: don't filter

% debug options
% set defaults first
constants.numSources = 5;
dbDeltaInitPos = [];
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
        % set perturbation of starting positions -> better perturb outside
        if isfield(dbOpt,'deltaInitPos')
            dbDeltaInitPos = dbOpt.deltaInitPos;
        end
        % set gradient option -> move into dataProperties
        if isfield(dbOpt,'gradientOption')
            constants.gradientOption = dbOpt.gradientOption;
        end
        % data for the analysis of F-statistics:
        % structure array (nSources x nTags x nTargets)
        % .residuals  vector of residuals
        % .resImSize  size of residual image
        % .initResiduals  vector of residuals of initial guess
        % .initResImSize  size of initial residual image
        % .source     source image
        % .sourceImSize size of source image
        % .sigma0     ssq of residuals
        % .fProb      p-value of f-test
        % .noise      estimated noise
        % .df         [n1, n2] degrees of freedom
        % .info   ct, sourceT, targetT, exitflag, success
        if isfield(dbOpt,'fStats')
            debugData.fStats = [];
        end
        % tbd
        if isfield(dbOpt,'objectiveFunction')
            debugData.objectiveFunction = 1;
        end
        % track results:
        % structure array (nSources x nTags x nTargets)
        % .startEndDelta 2x3 array with start delta, end delta
        % .sigma0 2x1 vector with start/end ssq of residuals
        % .info   ct, sourceT, targetT, exitflag, success
        % (in caller: add fields like deltaStartEndDelta as comparison to
        % reality)
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
    % searchRadius: sqrt of Q * factor or from idlist.trackInit
    searchRadius = sqrt(reshape(squeeze(qMatDiag)',[3,nTags])') *...
        constants.trackerRadiusMultiplicator;
    if ~isempty(idlist(t).trackInit)
        searchRadius(idlist(t).trackInit(:,1),:) = ...
            idlist(t).trackInit(:,2:end);
    end
    searchRadius4Source(t,:) = reshape(searchRadius',[1,3*nTags]);


    % for estimated tags: use trackInit


    % Fill in coords for all frames: We'll need them as starting
    % positions or for search radii.

    % coords: convert to pixels and switch x,y
    coords = idlist(t).linklist(:,[10,9,11])./p2m; %nTags x 3
    inputCoords(t,:,:) = ...
        reshape(coords,[1,nTags,3]); %[1, nTags, 3]
    coords4Source(t,:) = reshape(coords',[1,3*nTags]);
    inputAmp(t,:,:) = reshape(idlist(t).linklist(:,8),[1,nTags,1]);

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

% fill in source data
for iSource = sourceList'
    for iTag = 1:nTags
        fittingMatrices(iTag).A(iSource,iSource,1) = 1;
        fittingMatrices(iTag).B(iSource,1,:) = inputCoords(iSource,iTag,:);
        fittingMatrices(iTag).V(iSource,1,:) = inputQmatrixDiag(iSource,iTag,:);
    end
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
for t = sourceList'

    if iscell(movie)
        movieFrame = cdLoadMovie(movie,'',t);
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

end

%================================


%================================
%% TRACKING LOOP
%================================

% loop through the source/target pairs in trackingStrategy. For every
% target, loop until we have enough sources

% fittingIdx indicates the line where we will continue to write the fitting
% matrices
fittingIdx = sourceList(end);


% DEBUG
if debug && isfield(debugData,'trackResults')

    % track results:
    % structure array (nSources x nTags x nTargets)
    % .startEndDelta 2x3 array with start delta, end delta
    % .sigma0 2x1 vector with start/end ssq of residuals
    % .info   ct, sourceT, targetT, exitflag, success
    % (in caller: add fields like deltaStartEndDelta as comparison to
    % reality)

    tmp(1:nSources,1:nTags,1:constants.numSources) = ...
        struct('startEndDelta',zeros(2,3),'sigma0',zeros(2,1),...
        'info',zeros(1,5));
    debugData.trackResults = tmp;
    clear tmp
end
if debug && isfield(debugData,'fStats')

    % data for the analysis of F-statistics:
    % structure array (nSources x nTags x nTargets)
    % .residuals  vector of residuals
    % .resImGoodIdx
    % .initResiduals  vector of residuals of initial guess
    % .initResGoodIdx
    % .source     source image
    % .sourceImSize size of source image
    % .statistics  fprob, sigma0, targetNoise, n1, n2
    % .initStatistics  fprob, sigma0, targetNoise, n1, n2
    % .info   ct, sourceT, targetT, exitflag, success

    tmp(1:nSources,1:nTags,1:constants.numSources) = ...
        struct('residuals',[],'resImGoodIdx',[],'initResiduals',[],...
        'initResImGoodIdx',[],'source',[],'sourceImSize',[],...
        'sigma0',[],'statistics',zeros(1,5),...
        'initStatistics',zeros(1,5),'info',zeros(1,5));
    debugData.fStats = tmp;
    clear tmp

    % fill source images
    [debugData.fStats(:,:,1).source] = ...
        deal(sourceInfo(sourceList,:).intList);
    [debugData.fStats(:,:,1).sourceImSize] = ...
        deal(sourceInfo(sourceList,:).size);

end


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
        movieFrame = cdLoadMovie(movie,'',currentTarget);
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
        % sourceIdx is needed for fStats and trackResults
        targetIdx = 0;
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
            % sourceIdx is needed for fStats and trackResults
            sourceIdx = currentSource == sourceList;
            targetIdx = targetIdx + 1;
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
                if debug
                    % store exitflag
                    exitFlag(currentTarget,iTry,1) = xfl;
                end

            case 'dom'
                % use simple optimization which calculates parms
                % individually for each tag

                % loop settings
                maxIter = 100; % maximum number of iterations
                maxDelta = 0.01; % if the tags move less than this, quit

                % prepare fitting loop
                parameters = initialParameters;
                oldParameters = zeros(size(initialParameters,1),maxIter);
                sigmaResidual2 = zeros(nTags,maxIter);
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
                    fProbInit(iTag,1)  = fcdf(fRatioInit(iTag,1),dofInit(iTag,1),movieFrameSize);

                    %                     %f-test
                    %                     fRatio = sigmaResidual2(iTag,iter+1)/movieNoise(currentTarget)^2;
                    %                     fProb  = fcdf(fRatio,dof,movieFrameSize/2);
                    %                     isSuccess(iTag) = fProb < 0.95;

                    % debug
                    if debug && isfield(debugData,'fStats')

                        % store init - stuff
                        debugData.fStats(sourceIdx,iTag,targetIdx).initResiduals = residuals;
                        debugData.fStats(sourceIdx,iTag,targetIdx).initResImGoodIdx = targetInfo(iTag).goodIdx;
                        debugData.fStats(sourceIdx,iTag,targetIdx).initStatistics = ...
                            [fProbInit(iTag), sigmaResidual2Init(iTag), movieNoise(currentTarget).^2,...
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
                    % calculate 3d distance
                    delta = squeeze(sqrt(sum(delta.^2,1)));
                    % check for all distances below threshold
                    convergeIdx = find(all(delta < maxDelta, 1));

                    % if there is a convergeIdx, write into exitFlag how
                    % far back we had to go to find the same set of
                    % parameters. Else, continue loop.
                    if ~isempty(convergeIdx)
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
        [dummy,targetInfo] = extractIntensities(sourceInfo(currentSource,:),...
            targetInfo,movieFrame,constants,parameters,constants.gradientOption);

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

                    debugData.fStats(sourceIdx,iTag,targetIdx).residuals = residuals;
                    debugData.fStats(sourceIdx,iTag,targetIdx).resImGoodIdx = targetInfo(iTag).goodIdx;
                    debugData.fStats(sourceIdx,iTag,targetIdx).statistics = ...
                        [fProb(iTag), sigmaResidual2(iTag), movieNoise(currentTarget).^2,...
                        dof(iTag),movieFrameSize];

                end


            end



        else
            isSuccess = 0;
            fProb = -1;
        end


        % debug
        if debug && isfield(debugData,'trackResults')
            for iTag = 1:nTags
                debugData.trackResults(sourceIdx,iTag,targetIdx).startEndDelta =...
                    [initialParameters((iTag-1)*3+1:iTag*3);...
                    parameters((iTag-1)*3+1:iTag*3)];
                debugData.trackResults(sourceIdx,iTag,targetIdx).sigma0(1) = ...
                    sigmaResidual2Init(iTag);

                if isSuccess
                    debugData.trackResults(sourceIdx,iTag,targetIdx).sigma0(2) = ...
                        sigmaResidual2(iTag);
                end
            end

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
                debugData.trackResults(sourceIdx,iTag,targetIdx).info = ...
                    [ct, currentSource,currentTarget, xfl, success];
            end
            if isfield(debugData,'fStats')
                debugData.fStats(sourceIdx,iTag,targetIdx).info = ...
                    [ct, currentSource,currentTarget, xfl, success];
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
            disp(sprintf('%s, %i, %1.4f\n',ans{isSuccess+1},xfl, max(fProb)))
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
                fittingMatrices(iTag).A(fittingIdx,currentTarget,1) = 1;
                fittingMatrices(iTag).A(fittingIdx,currentSource,1) = -1;
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

        end %if isSuccess

        % whether it was successful or not, we remove the current s/t pair
        % from the strategy

        currentStrategy(1,:) = [];

    end % loop currentStrategy


end % while ~isempty trackingStrategy

%====================================



%====================================
%% CALCULATE POSITIONS
%====================================

% at this point, we basically just have to run myLscov to get the
% xyz-coords of the tags in the target frames.
% before, we need to remove the empty cols in A, though. These will make up
% the list of frames that could not be tracked.

% since all A-matrices are the same, we can find the failIdx on a single
% one
failIdx = (all(~fittingMatrices(1).A(:,:,1),1));
goodIdx = (1:nTimepoints)';
goodIdx(failIdx) = [];

% because of the way the new lscov is coded, it cannot cope with zero-rows
% anymore -> remove them
failIdx = (all(~fittingMatrices(1).A(:,:,1),2));
goodRows = 1:size(fittingMatrices(1).A,1);
goodRows(failIdx) = [];

% prepare "output"
outputCoords = zeros(nTimepoints,nTags,3);
[outputQmatrixDiag,mse] = deal(zeros(nTimepoints,3,nTags));

% loop through dimensions and fit
for iTag = 1:nTags
    for dim=1:3
        [outputCoords(goodIdx,iTag,dim), outputQmatrixDiag(goodIdx,dim,iTag), mse(goodIdx,dim,iTag)] = ...
            myLscov(fittingMatrices(iTag).A(goodRows,goodIdx,1),...
            fittingMatrices(iTag).B(goodRows,1,dim),fittingMatrices(iTag).V(goodRows,1,dim));
    end
end

% !!! Adjust the variance matrix. For whatever reason, the input error
% estimate is better than the output. In other words, we claim that the
% sigma zero is better estimated by detector and tracker than by the
% fitting, so we divide the output variance by sigma zero (mse).
outputQmatrixDiag(goodIdx,:,:) = ...
    outputQmatrixDiag(goodIdx,:,:).^2./mse(goodIdx,:,:);

% figure,plot(outputQmatrixDiag);
% hold on, plot(squeeze(inputQmatrixDiag),'.')
% sa=squeeze(sum(abs(fittingMatrices.A),1));
% hold on, plot(outputQmatrixDiag.*sqrt(sa),'o' );
% for iTag = 1:nTags
%     sa=repmat((sum(abs(fittingMatrices(iTag).A),1))',[1,3]);
%     outputQmatrixDiag(:,:,iTag) = outputQmatrixDiag(:,:,iTag) .* (sa);
% end

% if debug: shorten the matrices so that it isn't insanely big
if debug && isfield(debugData,'fitStats')
    for iTag=1:nTags
        fittingMatrices(iTag).A = ...
            sparse(fittingMatrices(iTag).A(goodRows,:,1));
        fittingMatrices(iTag).B = ...
            fittingMatrices(iTag).B(goodRows,:,:);
        fittingMatrices(iTag).V = ...
            fittingMatrices(iTag).V(goodRows,:,:);
    end
    % all A-matrices are identical
    for iTag=2:nTags
        fittingMatrices(iTag).A = [];
    end
    debugData.fitStats.fittingMatrices = fittingMatrices;
    debugData.fitStats.mse = mse;
    debugData.fitStats.outputQmatrixDiags = outputQmatrixDiags;
    debugData.fitStats.goodIdx = goodIdx;
end

%===================================

%===================================
%% WRITE OUT IDLISTTRACK
%===================================

% loop through idlist. Fill in the new coords and q-matrices. Keep the old
% entry wherever we had a lasting problem

idlisttrack = idlist;

for t = [goodTimes';NaN,goodTimes(1:end-1)']
    if any(t(1)==goodIdx)

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
        [dummy, sourceList]=find(fittingMatrices(1).A(targetIdx,:,1)==-1);
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

x=([-filterSizeXZ(1)/2:filterSizeXZ(1)/2])./halfSigmaXZ(1);
z=([-filterSizeXZ(2)/2:filterSizeXZ(2)/2])./halfSigmaXZ(2);
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



