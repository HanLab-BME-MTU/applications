function [runInfo]=polyDepolyMap(imDir,anDir,nPairs2analyze,doNorm,winSize)
%POLYDEPOLYMAP creates a map of polymer kinetics using mass conservation
%
% SYNOPSIS: [runInfo]=polyDepolyMap(imDir,anDir,nPairs2analyze,doNorm,winSize)
%
% INPUT: imDir: image directory (if [], will query user).
%        anDir: analysis (project) directory (if [], will query user).
%        nPairs2analyse: number of image pairs to analyze (if [], will
%                        analyze all).
%        n2avg: number of pairs to use for temporal averaging (if [], no
%               temporal averaging will be done).
%        gamma: gamma correction (default is 1, or no correction). gamma<1
%               boosts low values; 0>gamma>1 suppresses low values.
%        doNorm: 1 if you want to do (or re-do) image normalization;
%                otherwise this time-consuming step can be skipped
%        winSize: length of measurement window (default is 17 pixels)
%
% OUTPUT: runInfo: structure containing directory information and various
%                  calculated parameters.
%         /analysis/turn/turnRaw         : four maps are stored per frame pair,
%                                          corresponding to the four
%                                          spatially distinct grids used
%                                          for the turnover calculation.
%         /analysis/turnSpAvgDir/mapMats : one map per frame pair,
%                                          calculated by averaging the four
%                                          maps in /turnRaw
%
% NOTES: anDirectory should be set up like a normal FSM project with /corr
% and /edge.  However, /corr can contain the usual format with flowTrack
% files or Cyrus's specialized format.  In this format, the loaded
% flowTrack file contains a structure called flowHistory, which contains
% all the tracks.  In either case, a flowTrack file is needed for every
% frame pair (i.e. Time Step Size = 1 in imKymoAnalysis)
%
% Edge tracking should be done prior to running.  If image normalization
% has not been run before (ie there is no /norm directory), it will be done
% here by normalizing all the frames in imDir between 0 and 1.
%
% USERNAME: kathomps
% DATE: 12-Feb-2007
%
%try
tic

global DEBUG__




% these will keep track of min and max values over the whole movie
movieMin=0; movieMax=0;

% check input of image directory; query user if []
if nargin<1 || isempty(imDir)
    runInfo.imDir=uigetdir(pwd,'Please select image directory');
    imDir=runInfo.imDir;
else
    runInfo.imDir=imDir;
    if ~ischar(runInfo.imDir) || ~isdir(runInfo.imDir)
        error('POLYDEPOLYMAP: imDir must be a valid directory name')
    end
end

% check input of analysis directory; query user if []
if nargin<2 || isempty(anDir)
    runInfo.anDir=uigetdir(runInfo.imDir ,'Please select project analysis directory');
    anDir=runInfo.anDir;
else
    runInfo.anDir=anDir;
    if ~ischar(runInfo.anDir) || ~isdir(runInfo.anDir)
        error('POLYDEPOLYMAP: imDir must be a valid directory name')
    end
end

% check input of nPairs2analyze
if nargin<3 || isempty(nPairs2analyze)
    runInfo.nPairs2analyze=[]; %default - do them all
end

% check input of doNorm
if nargin<4 || isempty(doNorm) || doNorm~=1
    doNorm=0;
end

% check input of winSize
if nargin<5 || isempty(winSize)
    runInfo.winSize=17; % default - empirical estimate
else
    runInfo.winSize=winSize;
end

% check whether cell masks exist.  they are necessary for calculating mean
% background intensity during normalization.  if the images have already
% been normalized, and somehow the masks do not exist, then the user should
% run edge tracking or create them using ones(lengthOfImage,widthOfImage).
cmDir=[runInfo.anDir filesep 'edge' filesep 'cell_mask'];
if ~isdir(cmDir)
    error('POLYDEPOLYMAP: cmDir does not exist')
else % if valid directory, check for tifs
    [listOfCellMasks] = searchFiles('.tif',[],cmDir,0);
    if isempty(listOfCellMasks) % if no tifs, give error
        error('POLYDEPOLYMAP: the cell mask directory is empty')
    end
end

% create corr directory if it does not exist
corrDir=[runInfo.anDir filesep 'corr'];
if ~isdir(corrDir)
    mkdir(corrDir);
end

% create norm directory if it does not exist
normDir=[runInfo.anDir filesep 'norm' filesep 'normMats'];
if ~isdir(normDir)
    mkdir(normDir);
end

% create directories for turnover data
runInfo.turnDir=[runInfo.anDir filesep 'turn_winSize_' num2str(winSize)];
if ~isdir(runInfo.turnDir)
    mkdir(runInfo.turnDir);
else
    rmdir(runInfo.turnDir,'s');
    mkdir(runInfo.turnDir);
end
% subdirectory for data before averaging (based on 4 grids per frame)
turnRawDir=[runInfo.turnDir filesep 'turnRaw'];
mkdir(turnRawDir);
% subdirectory for spatial average of the 4 grids per frame
turnSpAvgDir=[runInfo.turnDir filesep 'turnSpAvgDir'];
mkdir([turnSpAvgDir filesep 'mapMats']);


% ===EXTRACT STORED FLOW TRACK DATA========================
% check to see if corr contains flow tracks in either the imKymoAnalysis
% format or Cyrus's format (ie flowHistory structure)
% if no files exist, user should run imKymoAnalysis to track flow or copy over
% flowTrack in Cyrus format into corr directory

% searchFiles not case sensitive; if you change this, be careful that names conform to
% standards

% vecPxyVxyAllFms: nFrames x 1 cell containing the tail positions
% and x- and y-components of all vectors
listFilesLin = searchFiles('flowTrack(\w*)_(\w*).mat',[],corrDir,0);
listFilesCyrus = searchFiles('(\w*)FlowTrack.mat','_',corrDir,0);

if ~isempty(listFilesLin(:,1)) && isempty(listFilesCyrus(:,1)) % Lin format
    runInfo.nFlwTrcks=length(listFilesLin(:,1)); %number of flowtracks

    % fill cell with measured velocity info from all tracks
    vecPxyVxyAllFms=cell(runInfo.nFlwTrcks,1);
    for n=1:runInfo.nFlwTrcks
        fileName=[char(listFilesLin(n,2)) filesep char(listFilesLin(n,1))];
        flwTrck=load(fileName);
        vecPxyVxyAllFms{n}=[flwTrck.flowTrack.p{1,1} flwTrck.flowTrack.v{1,1}];
    end

elseif isempty(listFilesLin(:,1)) && ~isempty(listFilesCyrus(:,1)) % Cyrus format
    fileName=[char(listFilesCyrus(1,2)) filesep char(listFilesCyrus(1,1))];
    flwTrck=load(fileName);
    flwTrck=flwTrck.flowHistory;
    runInfo.nFlwTrcks=length(flwTrck); %number of flowtracks

    % fill cell with measured velocity info from all tracks
    vecPxyVxyAllFms=cell(runInfo.nFlwTrcks,1);
    for n=1:runInfo.nFlwTrcks
        vecPxyVxyAllFms{n}=[flwTrck(1,n).p flwTrck(1,n).v];
    end

else % corr is either empty, contains more than 1 Cyrus .mat file, or contains no Lin flowTrack files
    error('POLYDEPOLYMAP: flow tracks missing from /corr or have unknown format')
end

% get total number of images to avoid incorrect naming and create string
% for naming files with correct number of digits
[listOfImages] = searchFiles('.tif',[],runInfo.imDir);
runInfo.nImTotExist=length(listOfImages);
s=length(num2str(runInfo.nImTotExist));
strg=sprintf('%%.%dd',s);

% get image dimensions from the first image in the list
[runInfo.imL,runInfo.imW]=size(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));

% find the number of images to analyze
if isempty(runInfo.nPairs2analyze)
    runInfo.nIm2an=runInfo.nFlwTrcks+1;
else
    runInfo.nIm2an=min(runInfo.nPairs2analyze,runInfo.nFlwTrcks)+1;
end

% check that the number of images to analyze is less than or equal to the
% total number that exist
if runInfo.nImTotExist < runInfo.nIm2an
    error('POLYDEPOLYMAP: total number of images should be >= number of images used in analysis')
end

% ===NORMALIZE USING CELL MASKS TO 0-1=====================
% or skip this time-consuming step if files exist in normDir
listNormIm = searchFiles('norm_image',[],normDir);
if isempty(listNormIm) || doNorm==1
    normImgSeries(runInfo,runInfo.nImTotExist);
else
    disp('POLYDEPOLYMAP: images in /norm/normMats assumed to be normalized')
end

% get orMask, the mask of all area covered by cell during any frame
orMask=zeros(runInfo.imL,runInfo.imW);
for i=1:runInfo.nIm2an
    fileNameMask=[char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))];
    cmask=imread(fileNameMask);
    orMask=cmask|orMask;
end

% get roiMask from intersection of fieldGeom-polygon and orMask (fieldGeom generated by selecting a ROI in
% imKymoAnalysis), or if it doesn't exist, use orMask
fieldGeom=[corrDir filesep 'fieldGeom.mat'];
if exist(fieldGeom,'file')
    fieldGeom=load([corrDir filesep 'fieldGeom']);
    roiMask=zeros(runInfo.imL,runInfo.imW); % size of whole image
    roiMask=roipoly(roiMask,fieldGeom.fieldGeom{1,1}.bndX,fieldGeom.fieldGeom{1,1}.bndY) & orMask;
else % if no roi was pre-selected, use whole cell
    roiMask=orMask;
end


% get array containing all speeds from raw data over whole movie
% this info will be used to calculate the optimal window size
% it is also saved in turnDir if further analysis desired (e.g. histogram)
velocityList=zeros(100000,1);
counter=0;
for i=1:runInfo.nFlwTrcks
    % magnitude of all velocity vectors in frame i
    mag=sqrt(vecPxyVxyAllFms{i,1}(:,3).^2+vecPxyVxyAllFms{i,1}(:,4).^2);

    % get rid of any row corresponding to a NaN vector
    vecPxyVxyAllFms{i,1}(isnan(mag),:)=[];
    mag(isnan(mag))=[];

    % extend array if getting too big
    nVecThisFm=length(mag);
    if counter+nVecThisFm>10000
        velocityList=[velocityList; zeros(100000,1)];
    end

    % fill vector with vector magnitudes
    velocityList(counter+1:counter+nVecThisFm)=mag;
    counter=counter+nVecThisFm;
end
velocityList(counter+1:end)=[];
save([runInfo.turnDir filesep 'velocityList'],'velocityList');

% background region standard deviation is calculated during image normalization
% load bg region std from all the frames
% use the max std from all frames in window size calculation
if exist([anDir filesep 'norm' filesep 'bgStd.mat'])
    bgStd=load([anDir filesep 'norm' filesep 'bgStd.mat']);
    bgStd=max(bgStd.bgStd);
else
    error('POLYDEPOLYMAT: bgStd.mat is missing - run normImgSeries')
end

% ===CALL FUNCTION TO FIND OPTIMAL WINDOW SIZE HERE========
% function hasn't been finished yet; call will go here.  should somehow use
% mean(velocityList) and bgStd.  for now, use empirical estimate.
% =========================================================

% area of windows in first frame
A1=runInfo.winSize^2;

% get coordinates for rectangle circumscribing roi (integer # of window)
% top-left corner (TLC) and bottom-right corner (BRC) of the rectangle
[yTLC,xTLC,yBRC,xBRC]=findRoiRect(roiMask,runInfo.winSize);

% need to make three other overlapping grids by shifting half a window length
% col1: m x n windows
% col2: m-1 x n windows (shifted down)
% col3: m x n-1 windows (shifted right)
% col4: m-1 x n-1 windows (shifted down and right)
shiftDR=ceil(runInfo.winSize/2); % amount to shift down or to the right
shiftUL=runInfo.winSize-shiftDR; % amount to shift up or to the left
yTLC=[yTLC yTLC+shiftDR yTLC         yTLC+shiftDR];
xTLC=[xTLC xTLC         xTLC+shiftDR xTLC+shiftDR];
yBRC=[yBRC yBRC-shiftUL yBRC         yBRC-shiftUL];
xBRC=[xBRC xBRC         xBRC-shiftUL xBRC-shiftUL];

recLength=yBRC-yTLC+1;
recWidth=xBRC-xTLC+1;
runInfo.nBoxesL=recLength/runInfo.winSize;
runInfo.nBoxesW=recWidth/runInfo.winSize;

% make sure findRoiRect picked a region with integer number of winSize x winSize
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
    testMask=testMask./max(testMask(:));
    indxStr1=sprintf(strg1,grid);
    imwrite(testMask,[runInfo.turnDir,filesep,'roi4calc_grid',indxStr1,'.tif']);
end

% allRecPixY/X are the yx-coordinates of the first grid (largest) over
% which the turnover map should be calculated
[allRecPixY,allRecPixX]=find(recMask(:,:,1));

% construct grid for vector interp - pyxNeeded will be used to get spline
% over the whole roi rectangle.  here we limit the number of grid points in
% each direction to keep the number of points to be interpolated somewhat
% small.  this of course can have the effect of making the vector map too
% rough if the roi is very large.
% nGrdPts=50;
% py=repmat(linspace(0.5,recLength(grid)+0.5,nGrdPts)',[1 nGrdPts])+yTLC(grid)-1;
% px=repmat(linspace(0.5,recWidth(grid)+0.5,nGrdPts),[nGrdPts 1])+xTLC(grid)-1;
% pyxNeeded=[py(:) px(:)]; %in cropped coordinates, covers all of cropped im

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

% created matrix of tiled window
%[boxNumTile]=tileboxnum(zeros(recLength(grid),recWidth(grid)),runInfo.winSize);

edgePix=cell(1,runInfo.nIm2an-1);
for i=1:runInfo.nIm2an-1 % loop thru all frame pairs
    %------------------------------------------------------------------
    % READ CELL MASKS AND LOAD NORMALIZED IMAGES FOR FRAMES i AND i+1

    if i==1
        fileNameMask1=[char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))];
        cellMask1=double(imread(fileNameMask1));

        indxStr1=sprintf(strg,i);
        im1=load([normDir filesep 'norm_image' indxStr1]);
        im1=im1.normImg;
    else
        cellMask1=cellMask2;
        im1=im2;
    end

    fileNameMask2=[char(listOfCellMasks(i+1,2)) filesep char(listOfCellMasks(i+1,1))];
    cellMask2=double(imread(fileNameMask2));

    indxStr2=sprintf(strg,i+1);
    im2=load([normDir filesep 'norm_image' indxStr2]);
    im2=im2.normImg;

    edgePix{i}=find(bwmorph(cellMask1,'remove'));

    %------------------------------------------------------------------
    % INTERPOLATE VECTOR FIELD FOR CURRENT FRAME AND GET SPLINE

    %dilate the cell mask by the window size to use for vector field interp.
    cellMaskDil=bwmorph(cellMask1,'dilate',runInfo.winSize+1); %ceil(maxThisFrm(i)));
    %check whether interp pixels are in or out of dilated cell
    inOrOut=cellMaskDil(pixInd);
    %Interp: [py0 px0 py1 px1], (all velocities initialized as zero)
    Interp=[pyxNeeded,pyxNeeded];

    % extract the vectors for frm1
    vecPxyVxy=vecPxyVxyAllFms{i,1};
    % separate into known vector tail positions and y and x components
    pyxVecKnown=vecPxyVxy(:,2:-1:1);
    vyxVecKnown=vecPxyVxy(:,4:-1:3);

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

    %     vy=Interp(:,3)-Interp(:,1);
    %     vx=Interp(:,4)-Interp(:,2);
    %     vyxVecNeeded=[vy vx];

    %     %plot stuff to check it


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
    pixInCellFm2=reshape(cellMask2(tAllPixInd),[recLength(1) recWidth(1)]);

    % normalize all of fm i+1 to (0,1), interpolate transformed pixels, and
    % correct their intensities by returning them to original range. this
    % is done b/c imInterp cannot work with images that are not in (0,1).
    im2norm01=(im2-min(im2(:)))/(max(im2(:))-min(im2(:)));
    im2interp01=reshape(imInterp(im2norm01, tAllPixYX),[recLength(1) recWidth(1)]);
    im2interpIOuterRec=(im2interp01*(max(im2(:))-min(im2(:))))+min(im2(:));

    %------------------------------------------------------------------
    % GET VALUES FOR THE FOUR GRIDS
    polyDepolyMap4Grids=zeros(runInfo.imL,runInfo.imW,4);
    for g=1:4
        % get pixel coordinates of the center of each window in grid g
        centersX=(ceil(runInfo.winSize/2):runInfo.winSize:recWidth(g)-ceil(runInfo.winSize/2)+1)+xTLC(g)-1;
        centersY=(ceil(runInfo.winSize/2):runInfo.winSize:recLength(g)-ceil(runInfo.winSize/2)+1)+yTLC(g)-1;
        [cCent rCent]=meshgrid(centersX,centersY);

        % get non-integer coordinates of corners of each window in grid g
        cornersX=(0.5:runInfo.winSize:recWidth(g)+0.5)+xTLC(g)-1;
        cornersY=(0.5:runInfo.winSize:recLength(g)+0.5)+yTLC(g)-1;
        [cCorn rCorn]=meshgrid(cornersX,cornersY);

        % transform box corners for grid g according to spline
        tBoxCornerGridYX=[fnval(spy,[rCorn(:)'; cCorn(:)'])' fnval(spx,[rCorn(:)'; cCorn(:)'])']...
            + [rCorn(:) cCorn(:)];

        % calc how many pixels per window are within the cell in fm i
        nCellPixInBox1=imResample(cellMask1(yTLC(g):yBRC(g),xTLC(g):xBRC(g)),[1 1],[runInfo.winSize runInfo.winSize]); % fm 1, grid g

        % calc how many pixels per window are within the cell in fm i+1
        nCellPixInBox2=imResample(pixInCellFm2((yTLC(g)-yTLC(1)+1):(yBRC(g)-yTLC(1)+1),...
            (xTLC(g)-xTLC(1)+1):(xBRC(g)-xTLC(1)+1)),[1 1],[runInfo.winSize runInfo.winSize]);


        % iob (in/out/border) keeps track of whether a window is
        %   0 - completely outside the cell in both frames
        %   1 - on the cell border in either frame
        %   2 - inside in both frames
        % initialize
        iob=zeros(size(nCellPixInBox1));
        % fill in where border=1
        iob((nCellPixInBox1>0) | (nCellPixInBox2>0))=1;
        % fill in where inside=2
        iob((nCellPixInBox1==A1) & (nCellPixInBox2==A1))=2;

        % denom keeps track of the area by which to divide I2-I1 to get a
        % per-pixel poly/depoly value.
        denom=nan*zeros(size(nCellPixInBox1));

        % where window is on the border in either frame (iob=1)
        % denom(iob==1)=abs(nCellPixInBox1(iob==1)-nCellPixInBox2(iob==1));
        %%%%%%%% CHECK THIS ASSUMPTION! %%%%%%% i think it is correct,
        % because it gives an average per-pixel value.
        % i tried the above value for denom and it seemed to boost the
        % edge values too high.
        denom(iob==1)=A1;

        % where window is inside the cell in both frames (iob=2)
        denom(iob==2)=A1;

        % get total intensity for each window in frame i
        I1=imResample(im1(yTLC(g):yBRC(g),xTLC(g):xBRC(g)),[1 1],[runInfo.winSize runInfo.winSize]);

        % find area of transformed window in frame i+1
        A2=grid2boxCorners(reshape(tBoxCornerGridYX(:,1),size(rCorn)),reshape(tBoxCornerGridYX(:,2),size(cCorn)))';
        % reshape A2 from vector form into matrix
        A2=reshape(A2,size(I1));

        % get total intensity for each window in frame i+1
        I2=imResample(im2interpIOuterRec((yTLC(g)-yTLC(1)+1):(yBRC(g)-yTLC(1)+1),...
            (xTLC(g)-xTLC(1)+1):(xBRC(g)-xTLC(1)+1)),[1 1],[runInfo.winSize runInfo.winSize]);

        % this is the beauty of the algorithm, i think.  in essence, to
        % calculate I2 we find the *estimated* average intensity per pixel in I2 using
        % image interpolation (above, just before for
        % loop).  we multiply this estimate by the actual area to get the
        % total I2 (total gray values in the window, not just per pixel).
        I2=(I2/A1).*A2;

        % turnover is simply the difference in intensity levels between the
        % window in frame i and its deformed image in i+1, divided by the
        % number of pixels used in the calculation (denom).
        rawTurnover=(I2-I1)./denom;

        % save the gridded values in the raw turnover data directory
        indxStr1=sprintf(strg,i);
        indxStr2=sprintf(strg,g);
        save([turnRawDir filesep 'rawTurnover' indxStr1 '_' indxStr2 '.mat'],...
            'rawTurnover','rCent','cCent','cellMask1');
    end

    if exist(DEBUG__) & i==1 %runInfo.nIm2an-1
        % plot to check during debugging
        bigGridMask=imread([runInfo.turnDir filesep 'roi4calc_grid1.tif']);
        figure(1); imshow(bigGridMask);
        
        hold on
        % plot window centers from frame i and show vector to where they go
        % in frame i+1 (blue)
        scatter(cCent(:),rCent(:),'bx') % window centers
        tCenter=[fnval(spy,[rCent(:)'; cCent(:)'])' fnval(spx,[rCent(:)'; cCent(:)'])']...
            + [rCent(:) cCent(:)];
        quiver(cCent(:),rCent(:),tCenter(:,2)-cCent(:),tCenter(:,1)-rCent(:),0,'b','lineWidth',2)
        scatter(tCenter(:,2),tCenter(:,1),'b.')
        
        % plot window corners from frame i and show vector to where they go
        % in frame i+1 (red)
        scatter(cCorn(:),rCorn(:),'r+') % window corners
        quiver(cCorn(:),rCorn(:),tBoxCornerGridYX(:,2)-cCorn(:),...
            tBoxCornerGridYX(:,1)-rCorn(:),0,'r','lineWidth',2)

        % plot raw vectors in green
        quiver(pyxVecKnown(:,2),pyxVecKnown(:,1),vyxVecKnown(:,2),vyxVecKnown(:,1),0,'g','lineWidth',2)

        %     quiver(pyxNeeded(:,2),pyxNeeded(:,1),vx,vy,0,'-r')
        %     hold on
        %     scatter(pyxStillNeeded(:,2),pyxStillNeeded(:,1),'+')
        
    end

end


% INTERPOLATE RAW DATA AND SAVE SPATIAL AVERAGE
[listOfRawTurnFiles] = searchFiles('rawTurnover',[],turnRawDir);
nPairs=size(listOfRawTurnFiles,1)/4;

s=length(num2str(nPairs));
strg=sprintf('%%.%dd',s);

counter=0;
for p=1:nPairs
    % use array to hold all four interpolated maps before averaging
    polyDepolyMap4Grids=zeros(runInfo.imL,runInfo.imW,4);

    for g=1:4
        counter=counter+1;
        fileName=[char(listOfRawTurnFiles(counter,2)) filesep char(listOfRawTurnFiles(counter,1))];
        turnFile=load(fileName);
        % interpolate using imDataMap
        polyDepolyMap4Grids(:,:,g)=imDataMap([runInfo.imL runInfo.imW],[turnFile.rCent(:) turnFile.cCent(:)],turnFile.rawTurnover(:),'grid',[runInfo.winSize runInfo.winSize]);
    end
    % spatially average over the four maps
    polyDepoly=nanmean(polyDepolyMap4Grids,3);
    
    % keep track of min/max for the movie
    movieMin=min(movieMin,nanmin(polyDepoly(:)));
    movieMax=max(movieMax,nanmax(polyDepoly(:)));
    
    % save the movie as a .mat file
    indxStr=sprintf(strg,p);
    save([turnSpAvgDir filesep 'mapMats' filesep 'polyDepoly' indxStr '.mat'],'polyDepoly');
end

runInfo.runTime=toc;
runInfo.movieMin=movieMin;
runInfo.movieMax=movieMax;

save([runInfo.turnDir filesep 'runInfo'],'runInfo');
save([runInfo.turnDir filesep 'edgePix'],'edgePix');

% catch
%     disp(['Something went wrong while analyzing ' runInfo.anDir])
%     [msgstr,msgid]=lasterr
%     return
% end

