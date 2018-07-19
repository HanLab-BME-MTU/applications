function [sourceInfo, targetInfo] = extractIntensities(sourceInfo, targetInfo, movieFrame, constants, parameters, doGradient)
%EXTRACTINTENSITIES is a subroutine to the tagTracker
% It extracts the intensities of the potentially overlapping tags and
% stores for each tag coordinates and intensities.
% source/target info are nSources x nTags structures with fields
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
%   .background     : background of image. Calculated as mean of image with
%                       masked intensity distributions
%   .backgroundVar  : variance of background. Used to correct for the noise
%                       reduction caused by interpolation
%   .interpCorr     : correction to the residuals to offset noise reduction
%                       caused by interpolation. Only works for linear
%                       interpolation. (calculation is hard for cubic and
%                       probably impossible for spline)
%
% (more to come)

% To check:
% - cut-off (currently at 1/50 of Gauss)
% - ratioList: use cut-off?
% - background: ok?
% - speed?
% - remove one-time use fields?
% - intensities: update? (i.e. calculate ratioList with intensities of
% target frame?)
% - targetInfo is not really needed as input argument! We could just copy
% the sourceInfo structure for preassignment




% check input

% if empty target - make source -- makeTarget = 0
% if none empty   - make all new target -- makeTarget = 1

if isempty(targetInfo)
    makeTarget = 0;
else
    makeTarget = 1;
end
if nargin < 6 || isempty(doGradient)
    doGradient = 0;
end

% interpolation scheme. Cell array with {scheme, correctOrNot}
if ~isfield(constants,'interpolation')
    constants.interpolation = {'*cubic',1};
end

% number of spots
nTags = length(sourceInfo);

% size of frame
frameSize = size(movieFrame);

% nCoords is for indexing into catenated coordList later
nCoords = ones(nTags+1,1);

switch makeTarget
    case 0
        % find size of source: Create 3D Gauss distribution around center
        % coordinate, then cut off wherever intensity is below the threshold,
        % and discard coordinates that are outside the movie frame
        % NEW: The list of pixels inside the threshold is still generated,
        % but the whole image is used

        % currently, use 1/50 of Gauss intensity (~3 sigma) as cut-off.

        % Gauss mask is the same for all Tags
        gaussSize = roundOddOrEven(6*constants.filterSigma,'odd','inf');
        gauss = GaussMask3D(constants.filterSigma, ...
            gaussSize);
        % cut-off: 1/50 of Gauss with max=1
        inGaussIdx = find(gauss > 1/50);

        % to prepare calculating list of coords: make grid, put 0 in center
        [gaussGridX,gaussGridY,gaussGridZ] = ...
            ndgrid(1:gaussSize(1),1:gaussSize(2),1:gaussSize(3));
        % gaussCoords are the coordinates of the block around the Gauss
        nGaussCoords = prod(gaussSize);
        gaussCoords = [gaussGridX(:),gaussGridY(:),gaussGridZ(:)] - ...
            repmat(ceil(gaussSize/2),[nGaussCoords,1]);

        % gridCoords are the coords of the gauss. We need them to eliminate
        % out-of-frame indices
        nInGaussCoords = length(inGaussIdx);
        [gridX, gridY, gridZ] = ind2sub(gaussSize,inGaussIdx);
        inGaussCoords = [gridX, gridY, gridZ];
        inGaussCoords = inGaussCoords...
            - repmat(ceil(gaussSize/2),[nInGaussCoords,1]);




        for iTag = 1:nTags

            % calculate list of coordinates and of the gauss coordinates
            coordList = gaussCoords + ...
                repmat(sourceInfo(iTag).centerCoord,[nGaussCoords,1]);
            inGaussList = inGaussCoords  + ...
                repmat(sourceInfo(iTag).centerCoord,[nInGaussCoords,1]);
            % find coords out of frame - out of frame is already if within the
            % border half pixel.
            goodCoordIdx = find((coordList(:,1)>1 & ...
                coordList(:,1)<frameSize(1)) &...
                (coordList(:,2)>1 & ...
                coordList(:,2)<frameSize(2)) &...
                (coordList(:,3)>1 & ...
                coordList(:,3)<frameSize(3)));
            goodInGaussIdx = find((inGaussList(:,1)>1 & ...
                inGaussList(:,1)<frameSize(1)) &...
                (inGaussList(:,2)>1 & ...
                inGaussList(:,2)<frameSize(2)) &...
                (inGaussList(:,3)>1 & ...
                inGaussList(:,3)<frameSize(3)));

            sourceInfo(iTag).coordList = coordList(goodCoordIdx,:);
            sourceInfo(iTag).goodIdx  = inGaussIdx(goodInGaussIdx);
            % store size of source after correction for out of frame. We
            % need to rebuild the deltaImage into an array for taking the
            % gradient.
            sourceInfo(iTag).size = ...
                round(max(sourceInfo(iTag).coordList,[],1) - ...
                min(sourceInfo(iTag).coordList,[],1) + [1,1,1]);
            nCoords(iTag+1) = length(goodCoordIdx);
            
            % if we had to make the sourceImage smaller at the lower left
            % corner ([1,1,1]), we need to update the index list.
            deltaSize = min(sourceInfo(iTag).coordList)-min(coordList);
            if ~all(deltaSize==0)
                deltaIdx = deltaSize(3)*gaussSize(1)*gaussSize(2) + ...
                    deltaSize(2)*gaussSize(1)*(gaussSize(3)-deltaSize(3)) +...
                    deltaSize(1)*(gaussSize(2)-deltaSize(2))*(gaussSize(3)-deltaSize(3));
                sourceInfo(iTag).goodIdx = sourceInfo(iTag).goodIdx - deltaIdx;
            end

        end

        clear gauss inGaussCoords inGaussList coordList nGoodGaussCoord gridX gridY gridZ

    case 1 % makeTarget == 1



        for iTag = 1:nTags

            %====================
            % here, we create the new image based on the new parameters.
            %  With the translation model, the new template is exactly
            % the source image. With the affine transform model, we have to
            % create a new template based on the transformation
            %
            % if affineTransform
            %       makeNewImage(source, parameters)
            %       readNewCoordList
            % else % just translation
            %       readCoordList from source and add delta
            % end

            deltaCoord = parameters(3*(iTag-1)+1:3*iTag)';
            % store deltaCoord
            targetInfo(iTag).deltaCoord = deltaCoord;
            coordList = sourceInfo(iTag).coordList + ...
                repmat(deltaCoord,...
                [size(sourceInfo(iTag).coordList,1),1]);

            % wrap around coordinates. Coordinates that have to be
            % shifted will show up in badIdx as either +1 or -1
            % for wraparound: don't add 1 to maximum, because we don't want
            % to change e.g. 0.5 into 10.5 if max(frameSize)=10
            [coordList,badIdx] = ...
                wraparound(coordList,[1,1,1;frameSize]);

            % find coords out of frame
            goodCoords = find(all(badIdx == 0,2));

            % goodCoords is what belongs REALLY to gauss (idx into
            % coordList)
            targetInfo(iTag).goodIdx = ...
                intersect(goodCoords,sourceInfo(iTag).goodIdx);

            % however, we keep the full coordinate list to extract
            % intensities.
            targetInfo(iTag).coordList = coordList;
            nCoords(iTag+1) = size(coordList,1);


        end

end % switch makeTarget
clear coordList goodCoords sx sy sz cl

%===========================
% BACKGROUND
%===========================
% calculate background if necessary. Mask Tags with NaNs, then nanmean.
% Store in info(1).background
switch makeTarget
    case 0 % we have to calc BG anyway

        sourceInfo = background(sourceInfo,frameSize,movieFrame);

    case 1
        % every time we re-calculate the background, as we're moving the
        % matrices around. Background is not going to change a lot, but
        % doing it some other way is going to be difficult

        targetInfo = background(targetInfo,frameSize,movieFrame);

end


%===============================
% RATIOLIST, INTLIST
%===============================

% calculate contribution of Tag intensity to pixels, interpolate
% intensities (in subfunction)
switch makeTarget
    case 0
        sourceInfo = gaussRatio2(sourceInfo,nCoords,constants,...
            movieFrame-sourceInfo(1).background);
    case 1
        [targetInfo.amp] = sourceInfo.amp;
        targetInfo = gaussRatio2(targetInfo,nCoords,constants,...
            movieFrame-targetInfo(1).background);
end

%=====================================
% CALCULATE CORRECTION FOR RESIDUALS
%=====================================

% Any interpolation scheme will result in a reduced noise between pixel
% centers. If the SNR is low, this results in a massive deformation of the
% objective function. We calculate the expected reduction to the residuals
% as a function of the subpixel-position of the tag and the variance of the
% background.
% Only calculate correction if there is a need
if constants.interpolation{2}
    switch makeTarget
        case 0
            sourceInfo = ...
                correctResiduals(sourceInfo,constants.interpolation{1});
        case 1
            targetInfo = ...
                correctResiduals(targetInfo,constants.interpolation{1});
    end
end


%=====================================
% TARGET: DELTAIMAGE, IMAGE GRADIENT
%=====================================
if makeTarget
    for iTag = 1:nTags
        % deltaInt: intensity difference between source and target image.
        % We subtract the full images, even if they have been wrapped
        % around
        % CORRECT FOR BLEACHING: adjust intensities of source according to
        % the estimated bleaching. Amplitudes are not the actual (noisy)
        % spot intensities from the dector, but the estimates from the
        % exponential bleaching curve.
        if constants.correctBleaching
        bleachingCorrection = targetInfo(iTag).amp/sourceInfo(iTag).amp;
        else
            bleachingCorrection = 1;
        end
        targetInfo(iTag).deltaInt = ...
            targetInfo(iTag).intList - ...
            sourceInfo(iTag).intList * bleachingCorrection;
        
        

        % ADD CORRECTION TO RESIDUALS - SEE CORRECTRESIDUALS FOR DETAILS
        % (preserve sign of deltaIntensity)
        if constants.interpolation{2}
            targetInfo(iTag).deltaInt = ...
                sign(targetInfo(iTag).deltaInt).*...
                sqrt(...
                targetInfo(iTag).deltaInt.^2 + ...
                sourceInfo(iTag).interpCorr * bleachingCorrection^2 ...
                + targetInfo(iTag).interpCorr);
        end

        % gradient: gradient of difference image - should not be affected
        % by the correction of the residuals, as we just added constants
        switch doGradient
            case 0
                % do nothing
            case 1
                % get filtered gradient
                % instead of first filtering the image and then taking the
                % gradient, I convolve the image with the gradient of the
                % Gaussian, using fastGauss3D for edge correction. To speed
                % up things, I use 3 1D kernels instead of 1 3D kernel
                                

                diffImage =...
                    - reshape(targetInfo(iTag).intList,sourceInfo(iTag).size);
                
%                 %% filter and gradient at the same time. For whatever
%                 %% reason the code only works fine if x and y are swapped.
%                 %% I have no clue why.
%                 dy = fastGauss3D(...
%                     fastGauss3D(...
%                     fastGauss3D(diffImage,[],[],1,constants.GaussGradXY),...
%                     [],[],1,constants.GaussXY'),...
%                     [],[],1,constants.GaussZ);
%                 dx = fastGauss3D(...
%                     fastGauss3D(...
%                     fastGauss3D(diffImage,[],[],1,constants.GaussXY),...
%                     [],[],1,constants.GaussGradXY'),...
%                     [],[],1,constants.GaussZ);
%                 dz = fastGauss3D(...
%                     fastGauss3D(...
%                     fastGauss3D(diffImage,[],[],1,constants.GaussXY),...
%                     [],[],1,constants.GaussXY'),...
%                     [],[],1,constants.GaussGradZ);
%                 %dx = fastGauss3D(diffImage,[],[],1,constants.GaussGradX);
%                 %dy = fastGauss3D(diffImage,[],[],1,constants.GaussGradY);
%                 %dz = fastGauss3D(diffImage,[],[],1,constants.GaussGradZ);
% 
%                 % since the blkdiag we have to make in the fitFcn is
%                 % annoying when reading out GaussPixels, we do that here
%                 % already
% 
%                 % collect gradient
%                 imGradient = [dx(:),dy(:),dz(:)];

                % filter diffImage
                diffImage = fastGauss3D(...
                    diffImage,constants.GaussSigma,constants.GaussSize,1);

                % gradient
                [dx,dy,dz] = deal(zeros(size(diffImage)));
                % central difference
                dx(2:end-1,:,:) = ...
                    (diffImage(3:end,:,:)-diffImage(1:end-2,:,:))/2;
                dy(:,2:end-1,:) = ...
                    (diffImage(:,3:end,:)-diffImage(:,1:end-2,:))/2;
                dz(:,:,2:end-1) = ...
                    (diffImage(:,:,3:end)-diffImage(:,:,1:end-2))/2;
                % forward difference for edges
                dx([1,end],:,:) = ...
                    diffImage([2,end],:,:) - diffImage([1,end-1],:,:);
                dy(:,[1,end],:) = ...
                    diffImage(:,[2,end],:) - diffImage(:,[1,end-1],:);
                dz(:,:,[1,end]) = ...
                    diffImage(:,:,[2,end]) - diffImage(:,:,[1,end-1]);
                
                imGradient = [dx(:),dy(:),dz(:)];

                % take only GaussPix
                if constants.gaussPixOnly
                    targetInfo(iTag).gradient = ...
                        imGradient(targetInfo(iTag).goodIdx,:);
                else
                    targetInfo(iTag).gradient = imGradient;
                end
                
            case 2
                % get gradient without filtering
                diffImage = ...
                    - reshape(targetInfo(iTag).deltaInt,sourceInfo(iTag).size);
                % gradient
                [dx,dy,dz] = deal(zeros(size(diffImage)));
                % central difference
                dx(2:end-1,:,:) = ...
                    (diffImage(3:end,:,:)-diffImage(1:end-2,:,:))/2;
                dy(:,2:end-1,:) = ...
                    (diffImage(:,3:end,:)-diffImage(:,1:end-2,:))/2;
                dz(:,:,2:end-1) = ...
                    (diffImage(:,:,3:end)-diffImage(:,:,1:end-2))/2;
                % forward difference for edges
                dx([1,end],:,:) = ...
                    diffImage([2,end],:,:) - diffImage([1,end-1],:,:);
                dy(:,[1,end],:) = ...
                    diffImage(:,[2,end],:) - diffImage(:,[1,end-1],:);
                dz(:,:,[1,end]) = ...
                    diffImage(:,:,[2,end]) - diffImage(:,:,[1,end-1]);
                
                % collect
                targetInfo(iTag).gradient = [dx(:),dy(:),dz(:)];
                
                % take only GaussPix if selected
                if constants.gaussPixOnly
                    targetInfo(iTag).gradient = ...
                        targetInfo(iTag).gradient(targetInfo(iTag).goodIdx,:);
                end
             
        end
    end
end

end % main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=========================================================================
% SUBFUNCTION BACKGROUND
%=========================================================================
function [infoStruct] = background(infoStruct,frameSize,bgImg)
% calculate background by masking Tag-pixels with NaNs, then nanmean

% pixel coords: 0.5-1.5 is within pixel 1
% since coordList is huge cubes, we should be fairly safe with respect to
% including unwanted intensity. Of course, there could always be some
% errant from another cell, but we assume that these stray signals are rare
% and their size small compared to the rest of the background pixels
allCoordPix = ceil(cat(1,infoStruct.coordList)-0.5);

% make index into image
allCoordIdx = sub2ind(frameSize,allCoordPix);

% make NaN-masked image (I believe this is faster than doing
% setdiff on the indices)
bgImg(allCoordIdx) = NaN;
bgImg = bgImg(:);

infoStruct(1).background = nanmean(bgImg);

% with the variance, I got an overestimation of the noise
% % also calculate and store variance of background signal for correcting the
% % noise-reduction due to interpolation
% % nanvar isn't in 6.51 - make compatible
% % infoStruct(1).backgroundVariance = nanvar(bgImg);
% infoStruct(1).backgroundVariance = var(bgImg(~isnan(bgImg)));

end % subfunction

%==========================================================================
% SUBFUNCTION GAUSSRATIO
%==========================================================================
% function [infoStruct] = gaussRatio(infoStruct,nCoords,constants,movieFrame)
%
% % to calculate the ratioMap, we will evaluate the Gauss of every Tag
% % at every gridpoint from which we'll later extract intensities.
%
% nTags = length(infoStruct);
% % catenate coordList
% allCoords = cat(1,infoStruct.coordList);
% % coords for the first Tag will be allCoords(allCoordsIdx(1):aci(2)-1)
% allCoordsIdx = [cumsum(nCoords)];
% allGauss = zeros(allCoordsIdx(end)-1, nTags);
%
% for iTag = 1:nTags
%
%     % get Gauss values for every center at every coordinate position
%     allGauss(:,iTag) = GaussListND(allCoords,...
%         constants.filterSigma,...
%         infoStruct(iTag).centerCoord) * infoStruct(iTag).amp;
%
% end
%
% % now calculate the ratioMap
% for iTag = 1:nTags
%
%     % get indexList into allGauss
%     idxList = [allCoordsIdx(iTag):allCoordsIdx(iTag+1)-1]';
%
%     % Gauss one col of allGauss, ratio is this col divided by the sum
%     % of all cols (for the idxList)
%     infoStruct(iTag).gaussList = allGauss(idxList,iTag);
%
%     % currently we divide by the sum of ALL the distributions. In the
%     % origninal code, there was no division below a certain threshold
%     infoStruct(iTag).ratioList = infoStruct(iTag).gaussList./...
%         sum(allGauss(idxList,:),2); % could put sum outside the loop
%
%     % with the ratioMap, we can go and read out intensities from the
%     % movie frames (which had their background removed in the function
%     % call)
%     %!!! interpolation uses the image coordinate system!!!
%     if ~strcmp(constants.interpolation{1},'*cubic')
%         infoStruct(iTag).intList = interp3(movieFrame,...
%             infoStruct(iTag).coordList(:,2),...
%             infoStruct(iTag).coordList(:,1),...
%             infoStruct(iTag).coordList(:,3),constants.interpolation{1}).*...
%             infoStruct(iTag).ratioList;
%     else
%         % use cubic interpolation directly
%         [nrows,ncols,npages] = size(movieFrame);
%         [s,t,w] = deal(zeros(size(infoStruct(iTag).coordList,1),1));
%         % bla
%         s(:) = infoStruct(iTag).coordList(:,2);
%         % bla
%         w(:) = infoStruct(iTag).coordList(:,3);
%         % bla
%         t(:) = infoStruct(iTag).coordList(:,1);
%
%         % Matrix element indexing
%         nw = (nrows+2)*(ncols+2);
%         ndx = floor(t)+floor(s-1)*(nrows+2)+floor(w-1)*nw;
%
%         % Compute intepolation parameters, check for boundary value.
%         d = find(s==ncols);
%         s(:) = (s - floor(s));
%         if length(d)>0, s(d) = s(d)+1; ndx(d) = ndx(d)-nrows-2; end
%
%         % Compute intepolation parameters, check for boundary value.
%         d = find(t==nrows);
%         t(:) = (t - floor(t));
%         if length(d)>0, t(d) = t(d)+1; ndx(d) = ndx(d)-1; end
%
%         % Compute intepolation parameters, check for boundary value.
%         d = find(w==npages);
%         w(:) = (w - floor(w));
%         if length(d)>0, w(d) = w(d)+1; ndx(d) = ndx(d)-nw; end
%
%         d = []; % Reclaim memory
%
%         % Expand v so interpolation is valid at the boundaries.
%
%             vv = zeros(size(movieFrame)+2);
%             vv(2:nrows+1,2:ncols+1,2:npages+1) = movieFrame;
%         vv(1,:,:)        = 3*vv(2,:,:)       -3*vv(3,:,:)     +vv(4,:,:); % Y edges
%         vv(nrows+2,:,:)  = 3*vv(nrows+1,:,:) -3*vv(nrows,:,:) +vv(nrows-1,:,:);
%         vv(:,1,:)        = 3*vv(:,2,:)       -3*vv(:,3,:)     +vv(:,4,:); % X edges
%         vv(:,ncols+2,:)  = 3*vv(:,ncols+1,:) -3*vv(:,ncols,:) +vv(:,ncols-1,:);
%         vv(:,:,1)        = 3*vv(:,:,2)       -3*vv(:,:,3)     +vv(:,:,4); % Z edges
%         vv(:,:,npages+2) = 3*vv(:,:,npages+1)-3*vv(:,:,npages)+vv(:,:,npages-1);
%         nrows = nrows+2; ncols = ncols+2; npages = npages+2;
%
%         % Now interpolate using computationally efficient algorithm.
%         F = zeros(size(s));
%         for iw = 0:3,
%
%             switch iw
%                 case 0
%                     ww = ((2-w).*w-1).*w;
%                 case 1
%                     ww = (3*w-5).*w.*w+2;
%                 case 2
%                     ww = ((4-3*w).*w+1).*w;
%                 case 3
%                     ww = (w-1).*w.*w;
%             end
%             for is = 0:3,
%                 switch is
%                     case 0
%                         ss = ((2-s).*s-1).*s;
%                     case 1
%                         ss = (3*s-5).*s.*s+2;
%                     case 2
%                         ss = ((4-3*s).*s+1).*s;
%                     case 3
%                         ss = (s-1).*s.*s;
%                 end
%                 for it = 0:3,
%                     switch it
%                         case 0
%                             tt = ((2-t).*t-1).*t;
%                         case 1
%                             tt = (3*t-5).*t.*t+2;
%                         case 2
%                             tt = ((4-3*t).*t+1).*t;
%                         case 3
%                             tt = (t-1).*t.*t;
%                     end
%                     F(:) = F + vv(ndx+(it+is*nrows+iw*nw)).*ss.*tt.*ww;
%                 end
%             end
%         end
%         F(:) = F/8;
%
%          infoStruct(iTag).intList = F(:).*...
%             infoStruct(iTag).ratioList;
%     end
%
% end
% end % subfunction

%==========================================================================
% SUBFUNCTION CORRECTRESIDUALS
%==========================================================================
function infoStruct = correctResiduals(infoStruct,interpolationString)

% to offset the reduction in residuals through the implicit filtering
% caused by interpolation, we calculate the expected reduction in residuals
% resulting from interpolation

% The residual image is the squared difference of two interpolated images.
% Therefore, it is not the expectation of the noise, but the squared
% expectation which equals the variance, that makes it into the residuals.
% The variance of interpolated noise with m=0 and var=V goes with the
% distance p from the center of one of two pixels for linear interpolation
% as: V(p)=(2p^2-2p+1)*V, becoming 1/2 the variance in the middle between
% two pixels (as to be expected for an average!). The difference to the
% variance at the pixels is 1-(...), and for 3D, the expression in the
% parentheses becomes the product of three independent functions of the
% subpixel shift in the respective dimensions.
%
% for cubic interpolation, V(p) = 0.25*(4-18p^2+16p^3+42p^4-60p^5+20p^6)
%
% A second source for (local) noise reduction is the splitting of
% intensities of overlapping tags - all the components to the signal should
% have the full noise added.

% find interpolation scheme
if strcmp(interpolationString(1),'*')
    interpolationScheme = interpolationString(2);
else
    interpolationScheme = interpolationString(1);
end

for iTag = 1:length(infoStruct)
    % get subpixel-shift
    subpixelShift = mod(infoStruct(iTag).coordList(1,:),1);
    % calculate product of reduction factor dependent on
    % interpolationScheme
    switch interpolationScheme
        case 'c' % cubic
            reductionFactor = prod(0.25*(...
                4-...
                18*subpixelShift.^2+...
                16*subpixelShift.^3+...
                42*subpixelShift.^4-...
                60*subpixelShift.^5+...
                20*subpixelShift.^6));

        case 'l' % linear
            reductionFactor = ...
                prod(2*subpixelShift.^2 - 2 * subpixelShift + 1);

    end
    % store interpCorr as the product of (1-reductionFactor) and the
    % variance with the noise reduction from the signal splitting factored
    % in
    infoStruct(iTag).interpCorr = ...
        ((1 - reductionFactor) * infoStruct(iTag).ratioList + ...
        (1 - infoStruct(iTag).ratioList)) * ...
        infoStruct(1).backgroundVariance;
end

end % subfunction









