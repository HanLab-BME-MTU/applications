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
[psfXY, psfZ] = calcFilterParms(...
    fitStruct.dataProperties.WVL,...
    fitStruct.dataProperties.NA,1.51,'bessel');
[gaussXY, gaussZ] = calcFilterParms(...
    fitStruct.dataProperties.WVL,...
    fitStruct.dataProperties.NA,1.51,'gauss');
add2psf = -0.1; % subtract 100 nm from grid, 0.3; % add 300 nm to grid
% npix_2 is always the same (grid goes from -n:n)
npix_2 = ceil((psfXY + add2psf/2)/pix2mu(1));


storeRatio = 0.1; % [0...1] ratio of acceptable NaN-pixels
nBackground = 3; % number of background pixels perpendicular to the axes

% initialize output. Careful: s2c2 is measured from c2!!!
intensities(1:nTimepoints) = struct('spotIntensities',NaN(4,1),...
    's1c1Int',[],...
    's1c1Vec',[],...
    's1c1VecPix',[],...
    's1c1Ang',[],...
    's1c1Bg' ,[],...
    's1c1Sp',[],...
    's2c2Vec',[],...
    's2c2VecPix',[],...
    's2c2Int',[],...
    's2c2Ang',[],...
    's2c2Bg' ,[],...
    's2c2Sp',[],...
    'c1c2Int',[],...
    'c1c2Vec',[],...
    'c1c2VecPix',[],...
    'c1c2Ang',[],...
    'c1c2Bg' ,[],...
    'c1c2Sp',[],...
    'nSpots',[],...
    's1s2Ang',[],...
    's1s2Vec',[]);

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
        rawMovie = cdLoadMovie({fitStruct.rawMovieName,'corr/raw'}, [], t);
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
        spotIdx(4) = find(strcmp(idlist(1).stats.labelcolor,'cen2'));
    end

    % update linklist, indices, according to order
    linklist = idlist(t).linklist(spotIdx,:);
    linklist(:,9:11) = linklist(:,[10,9,11]);

    % reset detectedTagIdx
    detectedTagIdx = find(linklist(:,2) > 0);
    nSpots = length(spotIdx);

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
        linklist(detectedTagIdx,8);
    intensities(t).nSpots = nSpots;




    % mask spots in rawMovie
    [maskMovie,bgMovie,wMovie] = deal(zeros(movieSize(1:3)));
    tmp = zeros([movieSize(1:3),nSpots]);

    % make a gaussMask for every single spot, then add up intensities and
    % backgrounds (it is conceivable that there are different backgrounds!)
    for i=1:nSpots %detectedTagIdx'
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

    for i = 1:3

        switch i
            case 1
                v_1 = linklist(2,9:11) - linklist(1,9:11);
                p_1 = linklist(1,9:11);
            case 2
                if nSpots == 3
                    v_1 = -linklist(2,9:11) + linklist(3,9:11);
                    p_1 = linklist(2,9:11);
                else
                    v_1 = -linklist(4,9:11) + linklist(3,9:11);
                    p_1 = linklist(4,9:11);
                end % check nSpots
                
            case 3
                if nSpots == 3
                    v_1 = [];
                else
                    v_1 = linklist(4,9:11) - linklist(2,9:11);
                end % check nSpots
                p_1 = linklist(2,9:11);
        end % switch i

        if ~isempty(v_1)

            % axis vector (don't scale yet)
            % number vectors 1,2,3 so that is easier to copy/paste

            [n_1,e_1] = normList(v_1);

            % angles
            angXY = atan2(e_1(2),e_1(1))*180/pi;
            angZ = 90 - acos(e_1(3))*180/pi;

            % xy-vector - normalize
            e_2 = [-e_1(2),e_1(1),0]/sqrt(sum(e_1(1:2).^2));

            % z-vector
            %e_3 = cross(e_1,e_2);
            e_3 = [0,0,1];

            % set vector length in microns. For axis, check how many pixels we get
            % in, approximately
            npix_1 = ceil(n_1/pix2mu(1));
            e_1 = e_1 * n_1/npix_1;
            e_2 = e_2 * pix2mu(1);
            e_3 = e_3 * pix2mu(1);

            % find how many pixels we need in perpendicular directions
            % grid goes from -n:n
            %npix_2 = ceil((psfXY + add2psf/2)/pix2mu(1));
            if ~all(e_3==[0,0,1])
                npix_3 = ceil((cos(angZ/180*pi)*psfZ+add2psf/2)/pix2mu(1));
            else
                npix_3 = ceil((psfZ + add2psf/2)/pix2mu(1));
            end

            % transform vectors to pixels
            e_1_pix = e_1./pix2mu;
            e_2_pix = e_2./pix2mu;
            e_3_pix = e_3./pix2mu;


            % make arbitraryGrid. Have npix_1 pixels along axis. This will
            % make the grid look like it's too small. However, we want to
            % look at intensities exactly from the spot positions on, and
            % not 0.5 pixels beyond.
            [xGrid,yGrid,zGrid] = arbitraryGrid(e_1_pix, e_2_pix, e_3_pix, ...
                p_1./pix2mu,...
                [0.5 npix_1-0.5], [-npix_2-nBackground npix_2+nBackground],...
                [-npix_3-nBackground npix_3+nBackground], movieSize(1:3));

            % make backgroundMask
            backgroundMask = false(size(xGrid));
            backgroundMask(:,[1:nBackground,end-nBackground+1:end],:) = true;
            backgroundMask(:,:,[1:nBackground,end-nBackground+1:end]) = true;

            % --- check grid with Imaris
            %             pd(3).XYZ = [xGrid(backgroundMask)*pix2mu(1),yGrid(backgroundMask)*pix2mu(2),zGrid(backgroundMask)*pix2mu(3)];
            %             pd(3).color = [0.7403, 0.7513, 0.7464, 0];pd(2).spotRadius =
            %             0.012;
            %             pd(2).XYZ = [xGrid(~backgroundMask)*pix2mu(1),yGrid(~backgroundMask)*pix2mu(2),zGrid(~backgroundMask)*pix2mu(3)];
            %             pd(2).color = [0.6403, 0.6513, 0.6464, 0];pd(2).spotRadius = 0.012;
            %             pd(1).XYZ = linklist(:,9:11);
            %             pd(1).color = [1,0,0,0];pd(1).spotRadius = 0.05;
            %             imarisPlot3(pd,pix2mu,cat(4,maskMovie,correctedImage,rawMovie-bgMovie))

            % read intensity
            interpImg = interp3(correctedImage,yGrid,xGrid,zGrid,'*linear',NaN);

            % split intensities into background and foreground
            % foreground: remove all background-zeros
            intensityImg = interpImg .* ~backgroundMask;
            intensityImg(:,all(backgroundMask(1,:,:),3),:) = [];
            intensityImg(:,:,all(backgroundMask(1,:,:),2)) = [];
            % background: Reassemble into 2-d array, where the 1st
            % dimension stays the axis vector
            backgroundImg = interpImg(backgroundMask);
            backgroundImg = reshape(backgroundImg,npix_1,[]);

            % also: read from maskMovie for kymograph
            mskImg = interp3(maskMovie,yGrid,xGrid,zGrid,'*linear',NaN);
            mskImg = mskImg .* ~backgroundMask;
            mskImg(:,all(backgroundMask(1,:,:),3),:) = [];
            mskImg(:,:,all(backgroundMask(1,:,:),2)) = [];
            
            rawImg = interp3(rawMovie,yGrid,xGrid,zGrid,'*linear',NaN);
            rawBg = rawImg(backgroundMask);
            rawBg = reshape(rawBg,npix_1,[]);
            rawImg = rawImg .* ~backgroundMask;
            rawImg(:,all(backgroundMask(1,:,:),3),:) = [];
            rawImg(:,:,all(backgroundMask(1,:,:),2)) = [];

            % check whether to store intensityImage
            store =  sum(~isfinite(intensityImg(:))) < storeRatio*numel(intensityImg);


            switch i
                case 1
                    % don't store if any NaN
                    if store
                        intensities(t).s1c1Int = intensityImg;
                        intensities(t).s1c1Bg = backgroundImg;
                        intensities(t).s1c1Sp = mskImg;
                        intensities(t).s1c1IntRaw = rawImg;
                        intensities(t).s1c1BgRaw = rawBg;
                    end
                    intensities(t).s1c1Vec = [e_1;e_2;e_3];
                    intensities(t).s1c1VecPix = [e_1_pix;e_2_pix;e_3_pix];
                    intensities(t).s1c1Ang = [angXY, angZ];
                case 2
                    % s2c2 is measured from c2!!
                    if store
                        intensities(t).s2c2Int = intensityImg;
                        intensities(t).s2c2Bg = backgroundImg;
                        intensities(t).s2c2Sp = mskImg;
                        intensities(t).s2c2IntRaw = rawImg;
                        intensities(t).s2c2BgRaw = rawBg;
                    end
                    intensities(t).s2c2Vec = [e_1;e_2;e_3];
                    intensities(t).s2c2VecPix = [e_1_pix;e_2_pix;e_3_pix];
                    intensities(t).s2c2Ang = [-angXY, -angZ];
                case 3
                    % don't store if any NaN
                    if store
                        intensities(t).c1c2Int = intensityImg;
                        intensities(t).c1c2Bg = backgroundImg;
                        intensities(t).c1c2Sp = mskImg;
                        intensities(t).c1c2IntRaw = rawImg;
                        intensities(t).c1c2BgRaw = rawBg;
                    end
                    intensities(t).c1c2Vec = [e_1;e_2;e_3];
                    intensities(t).c1c2VecPix = [e_1_pix;e_2_pix;e_3_pix];
                    intensities(t).c1c2Ang = [angXY, angZ];
            end

        end % if there is something to calculate
    end % loop 3 times
end % loop timepoints