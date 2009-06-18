function intensities = cdCheckStretching_readIntensities(fitStruct)
% readIntensities creates a grid between spb and cen tags and interpolates
% the intensities
%
% INPUT fitStruct: structure with fields
%		-idlist or slist: If idlist, code will generate slist from it
%		-dataProperties
%       -movieDir
%		-rawMovieName
%
% OUTPUT intensities: structure of length(idlist)
%       - spotIntensities [spb1 cen1 spb2 cen2]
%       - s1c1Int
%       - s1c1Vec [e_1;e_2;e_3] - normed to pixelsize, not to 1!
%       - s1c1VecPix - vectors for reading from image - not perpendicular!
%       - s1c1Ang [angle to plane, angle to x]
%       - same for s2c2, c1c2
%
%
% if there is any cen1*, don't read out at all

idlist = fitStruct.idlist;
nTimepoints = length(idlist);
pix2mu = [fitStruct.dataProperties.PIXELSIZE_XY,...
    fitStruct.dataProperties.PIXELSIZE_XY,...
    fitStruct.dataProperties.PIXELSIZE_Z];
% [psfXY, psfZ] = calcFilterParms(...
%     fitStruct.dataProperties.WVL,...
%     fitStruct.dataProperties.NA,1.51,'bessel',...
%     fitStruct.dataProperties.sigmaCorrection);
[gaussXY, gaussZ] = calcFilterParms(...
    fitStruct.dataProperties.WVL,...
    fitStruct.dataProperties.NA,1.51,'gauss',...
    fitStruct.dataProperties.sigmaCorrection);
% add2psf = -0.1; % subtract 100 nm from grid, 0.3; % add 300 nm to grid
% npix_2 is always the same (grid goes from -n:n)
%npix_2 = ceil((psfXY + add2psf/2)/pix2mu(1));


storeRatio = 0.25; % [0...1] ratio of acceptable NaN-pixels
nBackground = [3,1]; % number of background pixels perpendicular to the axes

% initialize output.
intensities(1:nTimepoints) = struct('spotIntensities',NaN(4,1),...
    'residualInt',[],...
    'residualBg',[],...
    'rawInt',[],...
    'rawBg',[],...
    'maskInt',[],...
    'nSpots',[],...
    'direction',[],...
    'dispResidual',[],...
    'dispMsk',[],...
    'dispRaw',[],...
    'dispX',[]);


% loop timepoints
goodTimes = catStruct(1,'idlist.linklist(1)');
for t=goodTimes'
    % check for good tags
    detectedTagIdx = find(idlist(t).linklist(:,2) > 0 & idlist(t).linklist(:,5) < 2);
    estimatedTagIdx = find(idlist(t).linklist(:,3) == 1 & idlist(t).linklist(:,5) < 2);
    goodTagIdx = sort([detectedTagIdx;estimatedTagIdx]);

    spotIdx = zeros(4,1);

    % load movie
    if ~isnumeric(fitStruct.rawMovieName)
        rawMovie = cdLoadMovie({fullfile(fitStruct.movieDir,fitStruct.rawMovieName),'corr/raw'}, [], t);
    else
        rawMovie = fitStruct.rawMovieName;
    end
    movieSize = size(rawMovie);

    % label tags: check labelcolor
    if any(strcmp(idlist(1).stats.labelcolor,'?'))
        % find via measuring distances
        dist = pdist(idlist(t).linklist(goodTagIdx,9:11));
        [maxDist,maxDistIdx] = max(dist);
        [row,col]=find(tril(ones(length(goodTagIdx)),-1));
        % spbIdx: largest distance between tags - these are spbs
        spbIdx = [row(maxDistIdx);col(maxDistIdx)];
        % cenIdx : all other indices
        cenIdx = setdiff(goodTagIdx,spbIdx);

        % spotIdx is: spb1, cen1, spb2, cen2
        % spb1-cen1 is smallest distance
        scDm = distMat2(idlist(t).linklist(spbIdx,9:11),idlist(t).linklist(cenIdx,9:11));
        [minRow,minCol] = find(scDm == min(scDm(:)));
        % minRow is spb1, minCol cen1
        spotIdx(1) = spbIdx(minRow);
        spbIdx(minRow) = [];
        spotIdx(2) = cenIdx(minCol);
        cenIdx(minCol) = [];
        spotIdx(3) = spbIdx;
        if isempty(cenIdx)
            spotIdx(4) = []; % if empty, spotIdx becomes shortened
        else
            spotIdx(4) = cenIdx;
        end
    else
        spotIdx(1) = find(strcmp(idlist(1).stats.labelcolor,'spb1'));
        spotIdx(2) = find(strcmp(idlist(1).stats.labelcolor,'cen1'));
        spotIdx(3) = find(strcmp(idlist(1).stats.labelcolor,'spb2'));
        if length(idlist(1).stats.labelcolor) > 3
            spotIdx(4) = find(strcmp(idlist(1).stats.labelcolor,'cen2'));
        else
            spotIdx(4) = [];
        end
    end

    % update linklist, indices, according to order
    linklist = idlist(t).linklist(spotIdx,:);
    linklist(:,9:11) = linklist(:,[10,9,11]);
    nSpots = length(spotIdx);
    % reset detectedTagIdx
    detectedTagIdx = find(linklist(:,2) > 0);

    % if we're missing any spb, skip everything. If we're missing one cen,
    % make the other into Cen1. If we're missing both CENs, skip everything
    badTagIdx = find(linklist(:,2) == 0 | linklist(:,3)==1 | linklist(:,5)==2);

    goodFrame = true;
    llIdx = 1:nSpots;
    if any(ismember([1 3],badTagIdx))
        goodFrame = false;
    elseif (nSpots==3 && any(badTagIdx == 2)) || (nSpots==4 && all(ismember([2,4],badTagIdx)))
        goodFrame = false;
    end
    if goodFrame && ~isempty(badTagIdx) 
        if badTagIdx == 2
        % it's a 4 spot movie, we're missing CEN1
        linklist = linklist([1,4,3],:);
        nSpots = 3;
        llIdx = [1,3,2];
        elseif badTagIdx == 4
            linklist = linklist([1,2,3],:);
            nSpots = 3;
            llIdx = 1:3;
        else
            error('unhandled tag disappearance')
        end
    end

    if goodFrame


        % find angle of s1s2
        v_s1s2 = linklist(3,9:11) - linklist(1,9:11);
        [n_s1s2,e_s1s2] = normList(v_s1s2);
        % angles
        angXY = atan2(e_s1s2(2),e_s1s2(1))*180/pi;
        angZ = 90 - acos(e_s1s2(3))*180/pi;
        intensities(t).s1s2Ang = [angXY,angZ];
        intensities(t).s1s2Vec = e_s1s2;

        % fill in intensities
        intensities(t).spotIntensities(detectedTagIdx) = ...
            linklist(llIdx,8);
        intensities(t).nSpots = nSpots;




        % mask spots in rawMovie
        [maskMovie,bgMovie,wMovie] = deal(zeros(movieSize(1:3)));
        tmp = zeros([movieSize(1:3),nSpots]);

        % make a gaussMask for every single spot, then add up intensities and
        % backgrounds (it is conceivable that there are different backgrounds!)
        % only consider actual tags!
        for i=1:nSpots
            maskMovieTmp = GaussMask3D([gaussXY/pix2mu(1),gaussXY/pix2mu(1),gaussZ/pix2mu(3)],...
                movieSize(1:3),linklist(i,9:11)./pix2mu,0,1,1);
            bgMovie = bgMovie + maskMovieTmp * linklist(i,16);
            maskMovie = maskMovie + maskMovieTmp * linklist(i,8);
            wMovie = wMovie + maskMovieTmp;
            tmp(:,:,:,i) = maskMovieTmp;
        end
        % weighted mean of backgrounds
        bgMovie = bgMovie./wMovie;
        bgMovie(isnan(bgMovie)) = nanmin(bgMovie(:));

        correctedImage = rawMovie - maskMovie - bgMovie;


        % make the grids and fill the fields
        % grids: pixelsize is pixel_XY in the perpendicular directions, close
        % to that in axial direction (see normList ###Vec)
        % In xy/z-direction, the extent of the grid is roughly 2xpsf+300nm
        % Background will be read from the four edges parallel to the main axis
        % (4*2x2(or 1x1)xnAxis pixels)

        % Unfortunately, interpolation will skew the noise in the image, which
        % we need in order to make a statistical estimate of whether pixels
        % contain significant intensity.
        % Therefore, we need to read actual pixel values.
        % 1) Read 2-3 vectors (spb-cen-(cen)-spb).
        % 2) Decide on whether to read in x or y direction, calculate number of
        %    pixels.
        % 3) For each pixel along the spb axis: Read +/- 3 sigma of pixels
        %    centered around the pixel along the spb-cen-(cen)-spb line. Start
        %    reading at -3 sigma along the spb axis.
        % 4) Dimension ordering: r=z, c=perp, e3=spb axis

        % there may be movies where the direction doesn't make sense. Too bad.

        % collect vectors

        % read vectors from idlist
        [n_s1c1,e_s1c1] = normList(linklist(2,9:11) - linklist(1,9:11));

        if nSpots == 3
            [n_c1c2,e_c1c2] = deal([]);
            [n_s2c2,e_s2c2] = normList(-linklist(2,9:11) + linklist(3,9:11));
        else
            [n_c1c2,e_c1c2] = normList(linklist(4,9:11) - linklist(2,9:11));
            [n_s2c2,e_s2c2] = normList(-linklist(4,9:11) + linklist(3,9:11));
        end

        % make vector list (unit vectors in pixels)
        % s1s2/s1c1/c1c2/c2s2/s1s2
        vectorList = zeros(5,3);
        vectorList(1,:) = e_s1s2;
        vectorList(2,:) = e_s1c1;
        vectorList(4,:) = -e_s2c2;
        vectorList(5,:) = e_s1s2;
        if nSpots == 3
            vectorList(3,:) = [];
        else
            vectorList(3,:) = e_c1c2;
        end
        vectorList = vectorList./repmat(pix2mu,size(vectorList,1),1);


        % figure out the good direction: xy-vector within 0...45, 135...180°
        % means we're going along x-axis
        if abs(angXY)<45 || abs(angXY)>135
            direction = 1;
            permuteOrder = [3,2,1];
        else
            direction = 2;
            permuteOrder = [3,1,2];
        end

        % it couls be that there is a step backwards somewhere. In that case,
        % try the other direction. If that doesn't work, return NaN
        directionOk = false;
        lastChance = false;

        while ~directionOk

            % count extents
            n3Sigma = ceil([gaussXY gaussZ]*3./pix2mu([1 3])/...
                abs(cos(atan2(e_s1s2(3-direction),e_s1s2(direction)))));
            % add 4 sigma on both sides of the SPB, too
            oneSigma = rayleighFromOri(e_s1s2, fitStruct.dataProperties.WVL,...
                fitStruct.dataProperties.NA,1.51,'gauss',...
                fitStruct.dataProperties.sigmaCorrection);
            % startPos is exact in pix, startPix is rounded
            startPos = (linklist(1,9:11) - e_s1s2*4*oneSigma)./pix2mu;
            startPix = round(startPos);
            % make sure we start inside the movie
            startPix = max(startPix,[1,1,1]);
            startPix = min(startPix,movieSize(1:3));
            endPos = (linklist(3,9:11) + e_s1s2*4*oneSigma)./pix2mu;
            endPix = round(endPos);
            endPix = max(endPix,[1,1,1]);
            endPix = min(endPix,movieSize(1:3));

            % count pixels along the direction
            nPixDir = abs(endPix(direction) - startPix(direction) + 1);

            % collect positions
            posList = zeros(6,3);
            posList(1,:) = startPix; % pix to avoid rounding issues
            posList(2,:) = linklist(1,9:11)./pix2mu;
            posList(3,:) = linklist(2,9:11)./pix2mu;
            posList(5,:) = linklist(3,9:11)./pix2mu;
            posList(6,:) = endPix; % pix to avoid rounding issues
            if nSpots == 3
                posList(4,:) = [];
            else
                posList(4,:) = linklist(4,9:11)./pix2mu;
            end

            % scale unit vectors so that they advance a pixel in The Direction each
            % step
            vectorList = vectorList./repmat(vectorList(:,direction),1,3);

            % preassign images - preassign them as images with background, which
            % will be removed later. Dimensions: z,perp,direction
            [residualIntBg,rawIntBg,maskIntBg] = deal(...
                NaN((n3Sigma(2)+nBackground(2))*2+1,...
                (n3Sigma(1)+nBackground(1))*2+1,nPixDir));
            [residualInt,rawInt,maskInt] = deal(...
                NaN(n3Sigma(2)*2+1,...
                n3Sigma(1)*2+1,nPixDir));

            % check whether the direction is ok, i.e. whether the dir-coords are
            % monotonously increasing
            deltaPos = diff(posList(:,direction));
            directionOk = all(sign(deltaPos./deltaPos(1))>0);
            
            % it is possible that all directions are negative.
            dirSign = sign(deltaPos(1));


            if directionOk
                % check whether any delta is smaller than 0.5 and make it 0.5 by
                % shifting positions a little. Otherwise, two tags lie on the same
                % pixel
                deltaPosSmallIdx = find(abs(deltaPos)<0.5);
                posList(deltaPosSmallIdx+1,direction) = ...
                    posList(deltaPosSmallIdx,direction)+0.5;

                % loop and read pixel values
                currentVectorIdx = 1;
                currentPosIdx = 1;
                currentPosPix = posList(currentPosIdx,:); %currentPosPix: in int pix
                currentPos = currentPosPix; % startPos is rounded

                % debug
                coordStore = {};

                for i=1:nPixDir
                    % direction determines how we're reading data
                    switch direction
                        case 1
                            c1 = currentPosPix(1);
                            c2 = max(currentPosPix(2)-n3Sigma(1)-nBackground(1),1):...
                                min(currentPosPix(2)+n3Sigma(1)+nBackground(1),movieSize(2));
                            c3 = max(currentPosPix(3)-n3Sigma(2)-nBackground(2),1):...
                                min(currentPosPix(3)+n3Sigma(2)+nBackground(2),movieSize(3));
                            d1 = c3-(currentPosPix(3)-n3Sigma(2)-nBackground(2))+1;
                            d2 = c2-(currentPosPix(2)-n3Sigma(1)-nBackground(1))+1;
                        case 2
                            c1 = max(currentPosPix(1)-n3Sigma(1)-nBackground(1),1):...
                                min(currentPosPix(1)+n3Sigma(1)+nBackground(1),movieSize(1));
                            c2 = currentPosPix(2);
                            c3 = max(currentPosPix(3)-n3Sigma(2)-nBackground(2),1):...
                                min(currentPosPix(3)+n3Sigma(2)+nBackground(2),movieSize(3));
                            d1 = c3-(currentPosPix(3)-n3Sigma(2)-nBackground(2))+1;
                            d2 = c1-(currentPosPix(2)-n3Sigma(1)-nBackground(1))+1;
                    end
                    
                    good1 = d1>=1&d1<=(n3Sigma(2)+nBackground(2))*2+1;
                    d1(~good1) = [];
                    good2 = d2>=1&d2<=(n3Sigma(1)+nBackground(1))*2+1;
                    d2(~good2) = [];

                    % read intensities
                    tmpInt = correctedImage(c1,c2,c3);
                    tmp = permute(tmpInt,permuteOrder);
                    residualIntBg(d1,d2,i) = tmp(good1,good2);
                    tmpInt = rawMovie(c1,c2,c3);
                    tmp = permute(tmpInt,permuteOrder);
                    rawIntBg(d1,d2,i) = tmp(good1,good2);
                    tmpInt = maskMovie(c1,c2,c3);
                    tmp = permute(tmpInt,permuteOrder);
                    maskIntBg(d1,d2,i) = tmp(good1,good2);

                    % update to next position
                    if i<nPixDir
                        nextPos = currentPos + vectorList(currentVectorIdx,:) * dirSign;
                    end

                    % check whether that would move us beyond a tag position
                    if i<nPixDir && round(dirSign*nextPos(direction))+0.5 > dirSign*posList(currentPosIdx+1,direction)
                        % update currentPosIdx, currentVectorIdx, recalculate nextPos
                        currentPosIdx = currentPosIdx + 1;
                        currentVectorIdx = currentVectorIdx + 1;
                        % nextPos is the rounded nextPos
                        nextPos = posList(currentPosIdx,:);% + ...
                        %                 vectorList(currentVectorIdx,:)*...
                        %                 (ceil(posList(currentPosIdx,direction)) - ...
                        %                 posList(currentPosIdx,direction));
                    end
                    currentPos = nextPos;
                    currentPosPix = round(currentPos);

                    % debug
                    %                 [xx,yy,zz] = ndgrid(c1,c2,c3);
                    %                 coordStore{end+1} = [xx(:),yy(:),zz(:)];

                end

                % if direction is not ok:
            else
                % check whether we already switched
                if lastChance
                    % just claim the direction is good - the loop won't repeat, and
                    % the arrays have already been initalized as NaN
                    directionOk = true;
                else
                    % switch direction
                    direction = 3-direction;
                    permuteOrder = [6 3 3] - permuteOrder;
                    % remember that this is the last chance
                    lastChance = true;
                end
            end % if directionOk

        end % while ~directionOk

        % debug
        allCoord = cat(1,coordStore{:});
        if ~isempty(allCoord) && ~all(isnan(residualIntBg(:)))
            figure('Name',sprintf('t: %i',t))
            for z=1:16
                subplot(4,4,z)
                imshow(rawMovie(:,:,z),[]);
                zIdxL = allCoord(:,3) == z;
                hold on
                plot(allCoord(zIdxL,2),allCoord(zIdxL,1),'.','MarkerSize',1)
                posIdxL = round(posList(:,3))==z;
                plot(posList(posIdxL,2),posList(posIdxL,1),'.r')
            end
        end

        % read and remove background
        backgroundIdxL = false(size(residualIntBg));
        backgroundIdxL([1:nBackground(2),end-nBackground(2)+1:end],:,:) = true;
        backgroundIdxL(:,[1:nBackground(1),end-nBackground(1)+1:end],:) = true;

        residualBg = residualIntBg(backgroundIdxL);
        residualInt(:) = residualIntBg(~backgroundIdxL);
        rawBg = rawIntBg(backgroundIdxL);
        rawInt(:) = rawIntBg(~backgroundIdxL);
        maskInt(:) = maskIntBg(~backgroundIdxL);

        % write coordinate data.
        intensities(t).spotPos = posList;
        intensities(t).direction = direction;

        % write intensities-structure if there is more than x% not-nan
        if sum(isnan(residualInt(:))) < storeRatio*numel(residualInt)
            intensities(t).residualInt = residualInt;
            intensities(t).residualBg = residualBg;
            intensities(t).rawInt = rawInt;
            intensities(t).rawBg = rawBg;
            intensities(t).maskInt = maskInt;




            % actually, for display purposes we need to interpolate the data.
            % Otherwise, sampling differs too much
            % also, only display where we measure!

            % make unit vectors. SPB-vector has 25nm-steps. Make in pix already
            unitStep = 0.025;
            zFactor = 3;
            e_1 = e_s1s2 * unitStep ./pix2mu; % since everything is in microns, we need to divide the unit by 40!
            e_2 = [-e_1(2) e_1(1), 0];
            e_2 = e_2/norm(e_2) *unitStep./pix2mu;
            e_3 = [0,0,unitStep*zFactor]./pix2mu; % z unit in 75 nm

            % check how many pixels we need in e_1: we want 500nm before and after
            n1=ceil(n_s1s2/unitStep);
            n1Ext = ceil(0.5/unitStep);

            % make sure we take enough pixels in directions 2 and 3: transform
            % cen-coords into new system, find max and min extents
            newCoord = linklist(:,9:11)-linklist(ones(nSpots,1),9:11); % subtract zero
            newCoord = newCoord./repmat(pix2mu,nSpots,1); % change to pix
            newCoord = [e_1(:),e_2(:),e_3(:)]\newCoord'; % transform
            newCoord = sign(newCoord).*ceil(abs(newCoord)); % round to +/- inf
            maxExtent = max(newCoord',[],1); % get max
            minExtent = min(newCoord',[],1); % get min

            % in 3-direction, we want to add 800 nm
            n3Ext = ceil(0.8/unitStep/zFactor);

            % make arbitrary grid.
            [xGrid,yGrid,zGrid] = arbitraryGrid(e_1, e_2, e_3, linklist(1,9:11)./pix2mu,...
                [-n1Ext n1+n1Ext], [minExtent(2)-n1Ext maxExtent(2)+n1Ext],...
                [minExtent(3)-n3Ext maxExtent(3)+n3Ext], movieSize(1:3));
%             [xx,dummy,dummy] = arbitraryGrid(e_1, e_2, e_3, linklist(1,9:11)./pix2mu,...
%                 [-n1Ext n1+n1Ext], [-0 0], [-0 0], movieSize(1:3));
            % interpolate. Subtract background already
            interpImg = interp3(correctedImage,yGrid,xGrid,zGrid,'*linear',NaN) - nanmedian(residualBg(:));
            mskImg = interp3(maskMovie,yGrid,xGrid,zGrid,'*linear',NaN);
            rawImg = interp3(rawMovie,yGrid,xGrid,zGrid,'*linear',NaN) - nanmedian(rawBg(:));

            % make projection already
            intensities(t).dispResidual = nansum(nansum(interpImg,2),3);
            intensities(t).dispMsk = nansum(nansum(mskImg,2),3);
            intensities(t).dispRaw = nansum(nansum(rawImg,2),3);
            intensities(t).dispX = (-n1Ext:n1+n1Ext)'*unitStep;%xx;
        end
    end % if goodFrame
end % loop time




%
%
%
%     for i = 1:3
%
%         switch i
%             case 1
%                 v_1 = linklist(2,9:11) - linklist(1,9:11);
%                 p_1 = linklist(1,9:11);
%             case 2
%                 if nSpots == 3
%                     v_1 = -linklist(2,9:11) + linklist(3,9:11);
%                     p_1 = linklist(2,9:11);
%                 else
%                     v_1 = -linklist(4,9:11) + linklist(3,9:11);
%                     p_1 = linklist(4,9:11);
%                 end % check nSpots
%
%             case 3
%                 if nSpots == 3
%                    v_1 = [];
%                 else
%                     v_1 = linklist(4,9:11) - linklist(2,9:11);
%                 end % check nSpots
%                 p_1 = linklist(2,9:11);
%         end % switch i
%
%         if ~isempty(v_1)
%
%             % axis vector (don't scale yet)
%             % number vectors 1,2,3 so that is easier to copy/paste
%
%             [n_1,e_1] = normList(v_1);
%
%             % angles
%             angXY = atan2(e_1(2),e_1(1))*180/pi;
%             angZ = 90 - acos(e_1(3))*180/pi;
%
%             % xy-vector - normalize
%             e_2 = [-e_1(2),e_1(1),0]/sqrt(sum(e_1(1:2).^2));
%
%             % z-vector
%             %e_3 = cross(e_1,e_2);
%             e_3 = [0,0,1];
%
%             % set vector length in microns. For axis, check how many pixels we get
%             % in, approximately. Then add 10 pixels on each side
%             addedPix = 20;
%             npix_1 = ceil(n_1/pix2mu(1));
%             e_1 = e_1 * n_1/npix_1;
%             e_2 = e_2 * pix2mu(1);
%             e_3 = e_3 * pix2mu(1);
%             npix_1 = npix_1 + addedPix;
%
%             % find how many pixels we need in perpendicular directions
%             % grid goes from -n:n
%             %npix_2 = ceil((psfXY + add2psf/2)/pix2mu(1));
%             if ~all(e_3==[0,0,1])
%                 npix_3 = ceil((cos(angZ/180*pi)*psfZ+add2psf/2)/pix2mu(1));
%             else
%                 npix_3 = ceil((psfZ + add2psf/2)/pix2mu(1));
%             end
%
%             % transform vectors to pixels
%             e_1_pix = e_1./pix2mu;
%             e_2_pix = e_2./pix2mu;
%             e_3_pix = e_3./pix2mu;
%
%
%             % make arbitraryGrid. Have npix_1 pixels along axis. This will
%             % make the grid look like it's too small. However, we want to
%             % look at intensities exactly from the spot positions on, and
%             % not 0.5 pixels beyond.
%             % correction 9/19/07: extend the grid to -4 pixels, so that we
%             % see all of the SPB (and so that we can judge the quality of
%             % spot masking)
%             % --- switch according to vector!! add pixels only on SPB side
%             [xGrid,yGrid,zGrid] = arbitraryGrid(e_1_pix, e_2_pix, e_3_pix, ...
%                 p_1./pix2mu,...
%                 [-addedPix/2+0.5 npix_1-addedPix/2-0.5], ...
%                 [-npix_2-nBackground npix_2+nBackground],...
%                 [-npix_3-nBackground npix_3+nBackground], movieSize(1:3));
%
%             % make backgroundMask
%             backgroundMask = false(size(xGrid));
%             backgroundMask(:,[1:nBackground,end-nBackground+1:end],:) = true;
%             backgroundMask(:,:,[1:nBackground,end-nBackground+1:end]) = true;
%
%             % --- check grid with Imaris
%             %             pd(3).XYZ = [xGrid(backgroundMask)*pix2mu(1),yGrid(backgroundMask)*pix2mu(2),zGrid(backgroundMask)*pix2mu(3)];
%             %             pd(3).color = [0.7403, 0.7513, 0.7464, 0];pd(2).spotRadius =
%             %             0.012;
%             %             pd(2).XYZ = [xGrid(~backgroundMask)*pix2mu(1),yGrid(~backgroundMask)*pix2mu(2),zGrid(~backgroundMask)*pix2mu(3)];
%             %             pd(2).color = [0.6403, 0.6513, 0.6464, 0];pd(2).spotRadius = 0.012;
%             %             pd(1).XYZ = linklist(:,9:11);
%             %             pd(1).color = [1,0,0,0];pd(1).spotRadius = 0.05;
%             %             imarisPlot3(pd,pix2mu,cat(4,maskMovie,correctedImage,rawMovie-bgMovie))
%
%             % read intensity
%             interpImg = interp3(correctedImage,yGrid,xGrid,zGrid,'*linear',NaN);
%
%             % split intensities into background and foreground
%             % foreground: remove all background-zeros
%             intensityImg = interpImg .* ~backgroundMask;
%             intensityImg(:,all(backgroundMask(1,:,:),3),:) = [];
%             intensityImg(:,:,all(backgroundMask(1,:,:),2)) = [];
%             % background: Reassemble into 2-d array, where the 1st
%             % dimension stays the axis vector
%             backgroundImg = interpImg(backgroundMask);
%             backgroundImg = reshape(backgroundImg,npix_1,[]);
%
%             % also: read from maskMovie for kymograph
%             mskImg = interp3(maskMovie,yGrid,xGrid,zGrid,'*linear',NaN);
%             mskImg = mskImg .* ~backgroundMask;
%             mskImg(:,all(backgroundMask(1,:,:),3),:) = [];
%             mskImg(:,:,all(backgroundMask(1,:,:),2)) = [];
%
%             rawImg = interp3(rawMovie,yGrid,xGrid,zGrid,'*linear',NaN);
%             rawBg = rawImg(backgroundMask);
%             rawBg = reshape(rawBg,npix_1,[]);
%             rawImg = rawImg .* ~backgroundMask;
%             rawImg(:,all(backgroundMask(1,:,:),3),:) = [];
%             rawImg(:,:,all(backgroundMask(1,:,:),2)) = [];
%
%             % check whether to store intensityImage
%             store =  sum(~isfinite(intensityImg(:))) < storeRatio*numel(intensityImg);
%
%
%             switch i
%                 case 1
%                     % don't store if any NaN
%                     if store
%                         intensities(t).s1c1Int = intensityImg;
%                         intensities(t).s1c1Bg = backgroundImg;
%                         intensities(t).s1c1Sp = mskImg;
%                         intensities(t).s1c1IntRaw = rawImg;
%                         intensities(t).s1c1BgRaw = rawBg;
%                     end
%                     intensities(t).s1c1Vec = [e_1;e_2;e_3];
%                     intensities(t).s1c1VecPix = [e_1_pix;e_2_pix;e_3_pix];
%                     intensities(t).s1c1Ang = [angXY, angZ];
%                 case 2
%                     % s2c2 is measured from c2!!
%                     if store
%                         intensities(t).s2c2Int = intensityImg;
%                         intensities(t).s2c2Bg = backgroundImg;
%                         intensities(t).s2c2Sp = mskImg;
%                         intensities(t).s2c2IntRaw = rawImg;
%                         intensities(t).s2c2BgRaw = rawBg;
%                     end
%                     intensities(t).s2c2Vec = [e_1;e_2;e_3];
%                     intensities(t).s2c2VecPix = [e_1_pix;e_2_pix;e_3_pix];
%                     intensities(t).s2c2Ang = [-angXY, -angZ];
%                 case 3
%                     % don't store if any NaN
%                     if store
%                         intensities(t).c1c2Int = intensityImg;
%                         intensities(t).c1c2Bg = backgroundImg;
%                         intensities(t).c1c2Sp = mskImg;
%                         intensities(t).c1c2IntRaw = rawImg;
%                         intensities(t).c1c2BgRaw = rawBg;
%                     end
%                     intensities(t).c1c2Vec = [e_1;e_2;e_3];
%                     intensities(t).c1c2VecPix = [e_1_pix;e_2_pix;e_3_pix];
%                     intensities(t).c1c2Ang = [angXY, angZ];
%             end
%
%         end % if there is something to calculate
%     end % loop 3 times
% end % loop timepoints