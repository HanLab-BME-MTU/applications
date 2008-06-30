function turnoverMap_CAW(runInfo,nImgs2Analyze,winL)


% imDir: image directory
% anDir: analysis directory
% nImgs2Analyze: number of image pairs to calculate (must be <= nImages-1)
% doNorm: 0 if image normalization has been done; 1 if it needs to be done
% winL: measurement window length (must be odd)
% methodStr: 'frame0grid' if poly events should be stored on same grid as depoly;
%            'frame1grid' if poly events should be stored on window center location in
%             second frame (see comment above where this flag is checked
%             for the rationale)




global DEBUG__
if isempty(DEBUG__)
    DEBUG__=0;
end

if nargin<1
    error('turnoverMap: Not enough input parameters')
end
if ~isstruct(runInfo)
    runInfo=struct;
end

if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('turnoverMap: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir] = formatPath(runInfo.anDir);
    [runInfo.imDir] = formatPath(runInfo.imDir);
end

% check input of nImgs2Analyze
if nargin<2 || isempty(nImgs2Analyze)
    runInfo.nImgs2Analyze=[]; % default - do them all
else
    runInfo.nImgs2Analyze=nImgs2Analyze;
end

% check input of winL
if nargin<3 || isempty(winL)
    runInfo.winL=11; % default - empirical estimate
else
    runInfo.winL=winL;
end


runInfo.methodStr='CAW'; % will plot on original grid


% normalized image directory
normDir=[runInfo.imDir filesep 'norm' filesep 'normMats'];
[listOfNormImgs]=searchFiles('norm_image',[],normDir,0);

% flow tracking directory, retrieve flow data
filtDir=[runInfo.anDir filesep 'corr' filesep 'filt'];
temp=load([filtDir filesep 'filteredFlow.mat']);
filteredFlow=temp.filteredFlow;


% vector coverage mask directory, get mask list
vecCovDir=[filtDir filesep 'vectorCoverageMask'];
[listOfVecCovMasks]=searchFiles('vectorCoverage',[],vecCovDir,0);

% create directories for turnover data
runInfo.turnDir=[runInfo.anDir filesep 'turn_' runInfo.methodStr filesep 'turn_winL_' num2str(runInfo.winL)];
if isdir(runInfo.turnDir)
    rmdir(runInfo.turnDir,'s');
end
turnRawDir=[runInfo.turnDir filesep 'raw'];
mkdir(turnRawDir);
mkdir([runInfo.turnDir filesep 'interp']);


% per-pixel variation due to noise is std of background divided by the
% square root of the number of pixels used to make the calculation. thus we
% can say that poly/depoly values that are above, say, 2x this value are
% significant.
temp=load([runInfo.imDir filesep 'norm' filesep 'normInfo.mat']);
normInfo=temp.normInfo;
runInfo.perPixNoiseVar=normInfo.bgStd/runInfo.winL;


% initialize fields to store min/max poly/depoly
runInfo.polyDepolyMovieMin=10^6;
runInfo.polyDepolyMovieMax=0;


% find the number of images to analyze
if isempty(runInfo.nImgs2Analyze)
    runInfo.nImgs2Analyze=runInfo.nImages;
else
    runInfo.nImgs2Analyze=nanmin(runInfo.nImgs2Analyze,runInfo.nImages);
end
s=length(num2str(runInfo.nImgs2Analyze));
strg=sprintf('%%.%dd',s);


% get roiMask from anDir (made useing polyDepolyChooseROI),or if it doesn't
% exist, use whole image
roiFile=[runInfo.anDir filesep 'polyDepolyROI.tif'];
if exist(roiFile,'file')
    roiMask=imread(roiFile);
else
    roiMask=ones(runInfo.imL,runInfo.imW);
end
roiMaskDil=bwmorph(roiMask,'dilate',runInfo.winL);

% pixels on cell boundary
temp=load([runInfo.anDir filesep 'edgePix.mat']);
edgePix=temp.edgePix;

%---------------------------------------

% area of windows in first frame
A1=runInfo.winL^2;

% get coordinates for rectangle circumscribing roi (integer # of window)
% top-left corner (TLC) and bottom-right corner (BRC) of the rectangle
[yTLC,xTLC,yBRC,xBRC,nWinL,nWinW]=findCalcRegion(roiMaskDil,runInfo.winL);


% need to make three other overlapping grids by shifting half a window length
% col1: m x n windows
% col2: m-1 x n windows (shifted down)
% col3: m x n-1 windows (shifted right)
% col4: m-1 x n-1 windows (shifted down and right)
shiftDR=ceil(runInfo.winL/2); % amount to shift down or to the right
shiftUL=runInfo.winL-shiftDR; % amount to shift up or to the left
yTLC=[yTLC yTLC+shiftDR yTLC         yTLC+shiftDR];
xTLC=[xTLC xTLC         xTLC+shiftDR xTLC+shiftDR];
yBRC=[yBRC yBRC-shiftUL yBRC         yBRC-shiftUL];
xBRC=[xBRC xBRC         xBRC-shiftUL xBRC-shiftUL];

recLength=yBRC-yTLC+1;
recWidth=xBRC-xTLC+1;
runInfo.nBoxesL=recLength/runInfo.winL;
runInfo.nBoxesW=recWidth/runInfo.winL;

% make sure findRoiRect picked a region with integer number of winL x winL
% windows
if mod(runInfo.nBoxesL,1)~=0 | mod(runInfo.nBoxesW,1)~=0
    error('POLYDEPOLYMAP: roiRect not made of integer number of windows');
end

% four turnover maps will be interpolated from data calculated on the four
% grids, and these will be averaged to get more accurate spatial average.
% here we write an image corresponding to the area around the ROI where
% the turnover map will be calculated relative to the cell.  check these if
% you suspect edge effects.
recMask=zeros(runInfo.imL,runInfo.imW,4);
strg1=sprintf('%%.%dd',length(num2str(4)));
for grid=1:4
    recMask(yTLC(grid):yBRC(grid),xTLC(grid):xBRC(grid),grid)=1;
    testMask=(recMask(:,:,grid)+roiMask);
    testMask=uint8(round(255.*testMask./max(testMask(:))));
    indxStr1=sprintf(strg1,grid);
    imwrite(testMask,[runInfo.turnDir,filesep,'roi4calc_grid',indxStr1,'.tif']);
end

% allRecPixY/X are the yx-coordinates of the first grid (largest) over
% which the turnover map should be calculated
[allRecPixY,allRecPixX]=find(recMask(:,:,1));


% construct grid for vector interp - pyxNeeded will be used to get spline
% over the whole roi rectangle. here we interpolate every 5 pixels.
nPixApart4Interp=5;
nGrdPtsY=ceil(recLength(1)/nPixApart4Interp);
nGrdPtsX=ceil(recWidth(1)/nPixApart4Interp);
py=repmat(linspace(0.5,recLength(1)+0.5,nGrdPtsY)',[1 nGrdPtsX])+yTLC(1)-1;
px=repmat(linspace(0.5,recWidth(1)+0.5,nGrdPtsX),[nGrdPtsY 1])+xTLC(1)-1;
pyxNeeded=[py(:) px(:)];

% get yx-coords of pixels in which the pts used for interp. fall
% pts on right and bottom edges of a pixel go in that pixel; pts on top and
% left end up in adjacent pixel (except for pts on the top or left edges of
% the image - those get put in adjacent image pixel)
indYX=ceil(pyxNeeded-0.5);
indYX(indYX==0)=1;
% convert to pixel indices
pixInd=xy2index(indYX(:,2),indYX(:,1),runInfo.imL,runInfo.imW,1);


%------------------------------------------------------------------
% ITERATE THROUGH FRAME PAIRS
for pairNumber=1:runInfo.nImgs2Analyze-1

    %------------------------------------------------------------------
    % READ CELL MASKS AND LOAD NORMALIZED IMAGES FOR FRAMES i AND i+1

    if pairNumber==1
        % mask0: dilated polyDepolyROI x vectorCoverageMask
        fileNameMask0=[char(listOfVecCovMasks(pairNumber,2)) filesep char(listOfVecCovMasks(pairNumber,1))];
        mask0=double(imread(fileNameMask0)).*roiMaskDil;

        % load first frame image
        fileNameImg0=[char(listOfNormImgs(pairNumber,2)) filesep char(listOfNormImgs(pairNumber,1))];
        img0=load(fileNameImg0);
        img0=img0.normImg;
    else
        mask0=mask1;
        img0=img1;
    end
    % mask1: dilated polyDepolyROI x vectorCoverageMask
    fileNameMask1=[char(listOfVecCovMasks(pairNumber+1,2)) filesep char(listOfVecCovMasks(pairNumber+1,1))];
    mask1=double(imread(fileNameMask1)).*roiMaskDil;

    % load second frame image
    fileNameImg1=[char(listOfNormImgs(pairNumber+1,2)) filesep char(listOfNormImgs(pairNumber+1,1))];
    img1=load(fileNameImg1);
    img1=img1.normImg;


    %------------------------------------------------------------------
    % INTERPOLATE VECTOR FIELD FOR CURRENT FRAME AND GET SPLINE

    %dilate the cell mask by the window size to use for vector field interp.
    cellMaskDil=bwmorph(mask0,'dilate',runInfo.winL+1);
    %check whether interp pixels are in or out of dilated cell
    inOrOut=cellMaskDil(pixInd);
    %Interp: [py0 px0 py1 px1], (all velocities initialized as zero)
    Interp=[pyxNeeded,pyxNeeded];

    % separate into known vector tail positions and y and x components
    pyxVecKnown=filteredFlow(1,pairNumber).pyxFilt;
    vyxVecKnown=filteredFlow(1,pairNumber).vyxFilt;

    % calculate the d0 and threshold (radius for sparse distance matrix)
    % for vectorFieldSparseInterp
    D=createDistanceMatrix(pyxVecKnown,pyxVecKnown);
    sD=sort(D,2); % col 1 all zeros b/c dist from point to itself, second col are NN distances
    d0=round(3*mean(sD(:,2))); %this is 3x the average NN distance
    thresh=3*d0;
    %interpolation step - only on points inside dilated cell
    Mi=vectorFieldSparseInterp([pyxVecKnown pyxVecKnown+vyxVecKnown],pyxNeeded(inOrOut,:),thresh,d0,[]);

    %after interp, we need to find values where we got NaNs
    %here we get the yx-coords for those points we still need
    pyxStillNeeded=Mi(isnan(Mi(:,3)),1:2);
    %now we get the rows in Mi where the NaNs are
    nanIndices=find(isnan(Mi(:,3)));
    %count the NaNs
    nNans=length(nanIndices);
    %make distance matrix from NaN points to known (raw) vectors
    distNans2Known=createDistanceMatrix(pyxStillNeeded,pyxVecKnown);
    %sort and keep the indices for the three closest points to each NaN pt
    [distNans2KnownSorted,IX]=sort(distNans2Known,2);
    IX=IX(:,1:5);
    %take average of the components of those three points
    vyxStillNeeded=[mean(reshape(vyxVecKnown(IX(:),1),nNans,[]),2) mean(reshape(vyxVecKnown(IX(:),2),nNans,[]),2)];
    Mi(nanIndices,3:4)=pyxStillNeeded+vyxStillNeeded;
    Interp(inOrOut,:)=Mi;

    %fit spline to interpolated vector field
    [spy,spx]=spFitVecField(Interp);

    %------------------------------------------------------------------
    % TRANSFORM EVERY PIXEL IN THE RECTANGLE CONTAINING THE LARGEST GRID
    % tAllPixYX is where the center of each pixel in the rect goes in fm 2,
    % measured in original image coordinates
    tAllPixYX=[fnval(spy,[allRecPixY'; allRecPixX'])' fnval(spx,[allRecPixY'; allRecPixX'])']...
        + [allRecPixY allRecPixX];
    tAllPixYX(tAllPixYX<1)=1;
    tAllPixYX(tAllPixYX(:,1)>runInfo.imL,1)=runInfo.imL;
    tAllPixYX(tAllPixYX(:,2)>runInfo.imW,2)=runInfo.imW;

    % get pixel indices for transformed pixels
    tAllPixInd=xy2index(ceil(tAllPixYX(:,2)),ceil(tAllPixYX(:,1)),runInfo.imL,runInfo.imW,1);
    % see which pixels are within the cell in fm i+1 using cellMask2
    pixInCellFm2=reshape(mask1(tAllPixInd),[recLength(1) recWidth(1)]);


    % arrange flow into input for vectorFieldSparseInterp
    PFvalues = [pyxVecKnown(:, 1), ...
        pyxVecKnown(:, 2), ...
        pyxVecKnown(:, 1) + vyxVecKnown(:, 1), ...
        pyxVecKnown(:, 2) + vyxVecKnown(:, 2)]; %this could be outside grid loop

    % generate smooth images for density sampling
    blurKernel = fspecial('gaussian', 11, 5); %truncated yes, but ...
    smoothImA = filterRegion(img0, mask0, blurKernel);
    smoothImB = filterRegion(img1, mask1, blurKernel);


    %------------------------------------------------------------------
    % GET VALUES FOR THE FOUR GRIDS
    polyDepolyMap4Grids=zeros(runInfo.imL,runInfo.imW,4);
    for g=1:4
        % get pixel coordinates of the center of each window in grid g
        centersX=(ceil(runInfo.winL/2):runInfo.winL:recWidth(g)-ceil(runInfo.winL/2)+1)+xTLC(g)-1;
        centersY=(ceil(runInfo.winL/2):runInfo.winL:recLength(g)-ceil(runInfo.winL/2)+1)+yTLC(g)-1;
        [cCent rCent]=meshgrid(centersX,centersY);

        % get non-integer coordinates of corners of each window in grid g
        cornersX=(0.5:runInfo.winL:recWidth(g)+0.5)+xTLC(g)-1;
        cornersY=(0.5:runInfo.winL:recLength(g)+0.5)+yTLC(g)-1;
        [cCorn rCorn]=meshgrid(cornersX,cornersY);


        % interpolate flow at corners
        vecF = vectorFieldSparseInterp(PFvalues, [rCorn(:), cCorn(:)], 22, 11, []);
        Py = zeros(size(rCorn));
        Py(:) = vecF(:, 1);
        Px = zeros(size(cCorn));
        Px(:) = vecF(:, 2);
        Vy = zeros(size(rCorn));
        Vy(:) = vecF(:, 3) - vecF(:, 1);
        Vx = zeros(size(cCorn));
        Vx(:) = vecF(:, 4) - vecF(:, 2);

        % set velocity at locations outside cell to zero (this is
        % intentional for the calculation)
        ind = sub2ind(size(mask0), clamp(ceil(Py(:)), 1, 512), ...
            clamp(ceil(Px(:)), 1, 512));
        outside = mask0(ind) ~= 1;
        Vy(outside) = 0;
        Vx(outside) = 0;

        % displace corners by flow and calculate areas
        rCornB = rCorn + Vy;
        cCornB = cCorn + Vx;
        A2 = quadrilateralAreas(rCornB, cCornB);

        % calculate densities for before (A) and after (B)
        densityA = sampleDensities(smoothImA, rCorn, cCorn);
        densityB = sampleDensities(smoothImB, rCornB, cCornB);

        % calculate masses and turnover
        % actually we'll normalize to box area (A1) so as to not scale with
        % window size
        massA = densityA;
        massB = densityB .* (A2 / A1);
        netTurnover = massB - massA;


        % save the gridded values in the raw turnover data directory
        indxStr1=sprintf(strg,pairNumber);
        indxStr2=sprintf(strg,g);
        save([turnRawDir filesep 'raw' indxStr1 '_' indxStr2 '.mat'],...
            'netTurnover','rCent','cCent');
    end
end

% INTERPOLATE RAW DATA AND SAVE SPATIAL AVERAGE
[listOfRawTurnFiles] = searchFiles('raw',[],turnRawDir,0);
nPairs=size(listOfRawTurnFiles,1)/4;
s=length(num2str(nPairs));
strg=sprintf('%%.%dd',s);

counter=0;
for pairNumber=1:nPairs
    % use array to hold all four interpolated maps before averaging
    polyDepolyMap4Grids=zeros(runInfo.imL,runInfo.imW,4);

    % mask0: dilated polyDepolyROI x vectorCoverageMask
    fileNameMask0=[char(listOfVecCovMasks(pairNumber,2)) filesep char(listOfVecCovMasks(pairNumber,1))];
    mask0=double(imread(fileNameMask0)).*roiMaskDil;

    % mask1: dilated polyDepolyROI x vectorCoverageMask
    fileNameMask1=[char(listOfVecCovMasks(pairNumber+1,2)) filesep char(listOfVecCovMasks(pairNumber+1,1))];
    mask1=double(imread(fileNameMask1)).*roiMaskDil;

    maskBoth=mask0 | mask1;

    for g=1:4
        counter=counter+1;
        fileName=[char(listOfRawTurnFiles(counter,2)) filesep char(listOfRawTurnFiles(counter,1))];
        turnFile=load(fileName);
        % interpolate using imDataMap
        polyDepolyMap4Grids(:,:,g)=imDataMap([runInfo.imL runInfo.imW],...
            [turnFile.rCent(:) turnFile.cCent(:)],turnFile.netTurnover(:),...
            'infLen',11,...
            'grid',[5 5],...
            'mask',maskBoth);
    end
    % spatially average over the four maps
    polyDepoly=nanmean(polyDepolyMap4Grids,3);

    runInfo.polyDepolyMovieMin=nanmin(runInfo.polyDepolyMovieMin,nanmin(polyDepoly(:)));
    runInfo.polyDepolyMovieMax=nanmax(runInfo.polyDepolyMovieMax,nanmax(polyDepoly(:)));

    % save the movie as a .mat file
    indxStr=sprintf(strg,pairNumber);
    save([runInfo.turnDir filesep 'interp' filesep 'polyDepoly' indxStr '.mat'],'polyDepoly');
end


save([runInfo.turnDir filesep 'runInfo'],'runInfo');


% -------------- subfunctions -------------------

function areas = quadrilateralAreas(rCorn, cCorn)

rP = rCorn(2:end,2:end) - rCorn(1:end-1,1:end-1);
cP = cCorn(2:end,2:end) - cCorn(1:end-1,1:end-1);
rQ = rCorn(2:end,1:end-1) - rCorn(1:end-1,2:end);
cQ = cCorn(2:end,1:end-1) - cCorn(1:end-1,2:end);

areas = 0.5 * crossProducts(rP, cP, rQ, cQ);
areas(nonConvex(rCorn, cCorn)) = NaN;



function concave = nonConvex(rCorn, cCorn)

rA = rCorn(1:end-1,2:end) - rCorn(1:end-1,1:end-1);
cA = cCorn(1:end-1,2:end) - cCorn(1:end-1,1:end-1);
rB = rCorn(2:end,2:end) - rCorn(1:end-1,2:end);
cB = cCorn(2:end,2:end) - cCorn(1:end-1,2:end);
rC = rCorn(2:end,1:end-1) - rCorn(2:end,2:end);
cC = cCorn(2:end,1:end-1) - cCorn(2:end,2:end);
rD = rCorn(1:end-1,1:end-1) - rCorn(2:end,1:end-1);
cD = cCorn(1:end-1,1:end-1) - cCorn(2:end,1:end-1);

concave = (crossProducts(rA, cA, rB, cB) < 0) | ...
    (crossProducts(rB, cB, rC, cC) < 0) | ...
    (crossProducts(rC, cC, rD, cD) < 0) | ...
    (crossProducts(rD, cD, rA, cA) < 0);



function xprod = crossProducts(rA, cA, rB, cB)

A = [cA(:)'; rA(:)'; zeros(1, numel(rA))];
B = [cB(:)'; rB(:)'; zeros(1, numel(rB))];
C = cross(A, B, 1);
xprod = zeros(size(rA));
xprod(:) = C(3,:);



function filteredIm = filterRegion(im, mask, kernel)

warning('off', 'MATLAB:divideByZero');

im(~mask) = 0;
filteredIm = imfilter(im, kernel);
W = imfilter(double(mask), kernel);
filteredIm = filteredIm ./ W;
filteredIm(~mask) = 0;



function density = sampleDensities(im, rCorn, cCorn)

cCorn(isnan(cCorn)) = -1;
rCorn(isnan(rCorn)) = -1;
dCorn = interp2(im, cCorn, rCorn, 'cubic');

corners = {dCorn(1:end-1,1:end-1), ...
    dCorn(1:end-1,2:end), ...
    dCorn(2:end,2:end), ...
    dCorn(2:end,1:end-1)};
sum = zeros(size(corners{1}));
N = zeros(size(corners{1}));

for m = 1:numel(corners)
    corn = corners{m};
    valid = ~isnan(corn);
    sum(valid) = sum(valid) + corn(valid);
    N(valid) = N(valid) + 1;
end

density = sum ./ N;



function clamped = clamp(unclamped, lowerLimit, upperLimit)

unclamped(unclamped < lowerLimit) = lowerLimit;
unclamped(unclamped > upperLimit) = upperLimit;
clamped = unclamped;
