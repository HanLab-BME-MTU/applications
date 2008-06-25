function [polyDepoly,accumY,accumX,runInfo]=turnoverMap(imDir,anDir,nPairs2analyze,doNorm,winL,polyLocFlag)

% threshold for vector outlier detection
runInfo.thresh=50;

% minimum number of neighbors for vector outlier detection
runInfo.nNeighborLowerLimit=3;

% directional deviation (degrees) allowed for a vector compared to its
% neighbors
runInfo.dirDev=15;

% how many points on each side of the measurement window to use
runInfo.nPtsPerSide=4; 

% imDir: image directory
% anDir: analysis directory
% nPairs2analyze: number of image pairs to calculate (must be <= nImages-1)
% doNorm: 0 if image normalization has been done; 1 if it needs to be done
% winL: measurement window length (must be odd)
% polyLocFlag: 0 if poly events should be stored on same grid as depoly; 1
%              if poly events should be stored on window center location in
%              second frame (see comment above where this flag is checked 
%              for the rationale)


%imDir='S:\scripps\analysis\thompson\Theriot\blebbistatinCells\050224D04Mf_231-266\images-renamed';
%anDir='S:\scripps\analysis\thompson\Theriot\blebbistatinCells\050224D04Mf_231-266\analysis';

%anDir='S:\scripps\analysis\thompson\Theriot\2008-06-12\050224D04Af\analysis';
%imDir='S:\scripps\analysis\thompson\Theriot\2008-06-12\050224D04Af\images';

%imDir=[pwd filesep 'exp1' filesep 'images' filesep 'maskedTifs'];
%anDir=[pwd filesep 'exp1' filesep 'analysis'];

global DEBUG__
if isempty(DEBUG__)
    DEBUG__=0;
end

%------------------------------------------------------------------
% CHECK USER INPUT & CREATE NEW DIRECTORIES

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
else
    runInfo.nPairs2analyze=nPairs2analyze;
end

% check input of doNorm
if nargin<4 || isempty(doNorm) || doNorm~=1
    doNorm=0; % doNorm must be 1 to run image normalization
end

% check input of winL
if nargin<5 || isempty(winL)
    runInfo.winL=11; % default - empirical estimate
else
    runInfo.winL=winL;
end

% check input of polyLocFlag
if nargin<6 || isempty(polyLocFlag)
    runInfo.polyLocFlag=0; % default - will plot on original grid
else
    runInfo.polyLocFlag=polyLocFlag;
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
runInfo.turnDir=[runInfo.anDir filesep 'turn_winL_' num2str(runInfo.winL)];
if isdir(runInfo.turnDir)
    rmdir(runInfo.turnDir,'s');
end
mkdir([runInfo.turnDir filesep 'raw']);
mkdir([runInfo.turnDir filesep 'interp']);
mkdir([runInfo.turnDir filesep 'vectorCoverageMask']);
if DEBUG__==1
    mkdir([runInfo.turnDir filesep 'overlay']);
end

% initialize fields to store min/max poly/depoly
runInfo.polyDepolyMovieMin=10^6;
runInfo.polyDepolyMovieMax=0;

%------------------------------------------------------------------
% EXTRACT STORED FLOW TRACK DATA

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
    disp('POLYDEPOLYMAP: flow tracks missing from /corr or have unknown format')
end

%------------------------------------------------------------------
% NORMALIZE IMAGES & GET ROI MASK

% get total number of images to avoid incorrect naming and create string
% for naming files with correct number of digits
[listOfImages] = searchFiles('.tif',[],runInfo.imDir);
runInfo.nImTotExist=length(listOfImages);
s=length(num2str(runInfo.nImTotExist));
strg=sprintf('%%.%dd',s);

% get image dimensions from the first image in the list
[runInfo.imL,runInfo.imW]=size(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
imL=runInfo.imL; imW=runInfo.imW;

% find the number of images to analyze
if isempty(runInfo.nPairs2analyze)
    runInfo.nIm2an=runInfo.nFlwTrcks+1;
else
    runInfo.nIm2an=nanmin(runInfo.nPairs2analyze,runInfo.nFlwTrcks)+1;
end

% check that the number of images to analyze is less than or equal to the
% total number that exist
if runInfo.nImTotExist < runInfo.nIm2an
    error('POLYDEPOLYMAP: total number of images should be >= number of images used in analysis')
end

% NORMALIZE images to 0-1, or skip this time-consuming step if files exist
% in normDir
listNormIm = searchFiles('norm_image',[],normDir);
if isempty(listNormIm) || doNorm==1
    [normInfo]=normImgSeries(runInfo,runInfo.nIm2an);
else
    disp('POLYDEPOLYMAP: images in /norm/normMats assumed to be normalized')
    temp=load([runInfo.anDir filesep 'norm' filesep 'normInfo.mat']);
    normInfo=temp.normInfo;
end

% per-pixel variation due to noise is std of background divided by the
% square root of the number of pixels used to make the calculation. thus we
% can say that poly/depoly values that are above, say, 2x this value are
% significant.
runInfo.perPixNoiseVar=normInfo.bgStd/runInfo.winL;

% get roiMask from fieldGeom-polygon (fieldGeom generated by selecting a ROI in
% imKymoAnalysis), or if it doesn't exist, use whole image
fieldGeom=[corrDir filesep 'fieldGeom.mat'];
if exist(fieldGeom,'file')
    fieldGeom=load([corrDir filesep 'fieldGeom']);
    roiMask=zeros(imL,imW); % size of whole image
    roiMask=roipoly(roiMask,fieldGeom.fieldGeom{1,1}.bndX,fieldGeom.fieldGeom{1,1}.bndY);
else % if no roi was pre-selected, use whole cell
    roiMask=ones(imL,imW);
end
runInfo.roiMask=roiMask;

% initialize cell to contain pixels on cell boundary
runInfo.edgePix=cell(1,runInfo.nIm2an-1);
filteredFlow(runInfo.nIm2an)=struct('pyxRaw',[],'vyxRaw',[],'pyxFilt',[],'vyxFilt',[]);
%------------------------------------------------------------------
% ITERATE THROUGH FRAME PAIRS
for pairNumber=1:runInfo.nIm2an-1
    
    indxStr1=sprintf(strg,pairNumber);   
    
    %------------------------------------------------------------------
    % READ CELL MASKS AND LOAD NORMALIZED IMAGES FOR FRAMES i AND i+1

    if pairNumber==1
        % load first frame mask
        fileNameMask0=[char(listOfCellMasks(pairNumber,2)) filesep char(listOfCellMasks(pairNumber,1))];
        cellMask0=double(imread(fileNameMask0)).*roiMask;
        % find cell edge pixels for first frame
        runInfo.edgePix{pairNumber}=find(bwmorph(cellMask0,'remove'));
        % load first frame image
        indxStr1=sprintf(strg,pairNumber);
        img0=load([normDir filesep 'norm_image' indxStr1]);
        img0=img0.normImg;
    else
        cellMask0=cellMask1;
        img0=img1;
    end
    % load second frame mask
    fileNameMask1=[char(listOfCellMasks(pairNumber+1,2)) filesep char(listOfCellMasks(pairNumber+1,1))];
    cellMask1=double(imread(fileNameMask1)).*roiMask;
    % find the edge pixels for second image based on mask
    runInfo.edgePix{pairNumber+1}=find(bwmorph(cellMask1,'remove'));
    % load second frame image
    indxStr2=sprintf(strg,pairNumber+1);
    img1=load([normDir filesep 'norm_image' indxStr2]);
    img1=img1.normImg;

    % smooth the images
    blurKernel=fspecial('gaussian', 11, 5);
    img0 = imfilter(img0, blurKernel,'symmetric');
    img1 = imfilter(img1, blurKernel,'symmetric');

    %------------------------------------------------------------------

    % extract the vectors for frm1, separate into known vector tail
    % positions and y and x components
    vecPxyVxy=vecPxyVxyAllFms{pairNumber,1};
    pyxVecKnown=vecPxyVxy(:,2:-1:1);
    vyxVecKnown=vecPxyVxy(:,4:-1:3);

    % remove vector outliers
    [pyxFiltered,vyxFiltered,vectorCoverageMask,runInfo]=removeVectorOutliers(runInfo,pyxVecKnown,vyxVecKnown);

    filteredFlow(pairNumber).pyxRaw=pyxVecKnown;
    filteredFlow(pairNumber).vyxRaw=vyxVecKnown;
    filteredFlow(pairNumber).pyxFilt=pyxFiltered;
    filteredFlow(pairNumber).vyxFilt=vyxFiltered;
    
    % get grid for calculating poly/depoly
    % findCalcRegion gives coordinates for a rectangle surrounding the
    % vectorCoverageMask, the sides of which are a multiple of gridSpace.
    % we need to add 1 to each side to allow for integer coordinates. for
    % example, if the rectangle were 15 pixels long, and gridSpace=5; we
    % want to end up with coordinates at 1,6,11,16.
    gridSpace=5;
    [yTLC,xTLC,yBRC,xBRC,nWinL,nWinW]=findCalcRegion(vectorCoverageMask,gridSpace);
    if yBRC<imL
        yBRC=yBRC+1;
    elseif yTLC>1
        yTLC=yTLC-1;
    else
        warning('TURNOVERMAP: border problem in grid to calculate poly/depoly');
    end

    if xBRC<imW
        xBRC=xBRC+1;
    elseif xTLC>1
        xTLC=xTLC-1;
    else
        warning('TURNOVERMAP: border problem in grid to calculate poly/depoly');
    end

    % these are the grid points for the poly/depoly calculation
    [cX0 cY0]=meshgrid(xTLC:gridSpace:xBRC,yTLC:gridSpace:yBRC);
    cY0=cY0(:); cX0=cX0(:);

    %------------------------------------------------------------------

    % here we interpolate the vector field on the grid points, as well as
    % the measurement window edge points

    d0=3*runInfo.meanNN; % d0 is weight for interpolation
    %thresh=3*d0; % threshold is radius for sparse distance matrix

    % interpolation step for CENTER points
    Mi=vectorFieldInterp([pyxFiltered pyxFiltered+vyxFiltered], [cY0, cX0], d0, []);
    %after interp, we need to find values where we got NaNs
    %here we get the rows in Mi where the NaNs are
    nanIndices=find(isnan(Mi(:,3)));
    %count the NaNs
    nNans=length(nanIndices);
    %here we get the yx-coords for those points we still need
    pyxStillNeeded=Mi(nanIndices,1:2);
    %make distance matrix from NaN points to known (raw) vectors
    distNans2Known=createDistanceMatrix(pyxStillNeeded,pyxFiltered);
    %sort and keep the indices for the three closest points to each NaN pt
    [distNans2KnownSorted,IX]=sort(distNans2Known,2);
    IX=IX(:,1:3);
    %take average of the components of those three points
    vyxStillNeeded=[mean(reshape(vyxFiltered(IX(:),1),nNans,[]),2) mean(reshape(vyxFiltered(IX(:),2),nNans,[]),2)];
    Mi(nanIndices,3:4)=pyxStillNeeded+vyxStillNeeded;
    % assign the transformed center point coordinates
    cY1=Mi(:,3);
    cX1=Mi(:,4);

    
    % get nEdgePts x nWindows arrays containing coordinates of the
    % points on the window perimeter (including corners)
    [edgeY0,edgeX0]=findEdgePts(cY0,cX0,runInfo.winL,runInfo.nPtsPerSide);
    [nEdgePts,nWindows]=size(edgeY0);
    counter=1;
    Mi=zeros(nEdgePts*nWindows,4);
    for iEdge=1:nEdgePts
        % interpolation step for window EDGE points
        Mi(counter:counter+nWindows-1,:)=vectorFieldInterp([pyxFiltered pyxFiltered+vyxFiltered], [edgeY0(counter:counter+nWindows-1)', edgeX0(counter:counter+nWindows-1)'], d0, []);
        counter=counter+nWindows;
    end
    %after interp, we need to find values where we got NaNs
    %here we get the rows in Mi where the NaNs are
    nanIndices=find(isnan(Mi(:,3)));
    %count the NaNs
    nNans=length(nanIndices);
    %here we get the yx-coords for those points we still need
    pyxStillNeeded=Mi(nanIndices,1:2);
    %make distance matrix from NaN points to known (raw) vectors
    distNans2Known=createDistanceMatrix(pyxStillNeeded,pyxFiltered);
    %sort and keep the indices for the three closest points to each NaN pt
    [distNans2KnownSorted,IX]=sort(distNans2Known,2);
    IX=IX(:,1:3);
    %take average of the components of those three points
    vyxStillNeeded=[mean(reshape(vyxFiltered(IX(:),1),nNans,[]),2) mean(reshape(vyxFiltered(IX(:),2),nNans,[]),2)];
    Mi(nanIndices,3:4)=pyxStillNeeded+vyxStillNeeded;
    % assign the transformed edge point coordinates
    edgeY1=reshape(Mi(:,3),4*runInfo.nPtsPerSide-3,[]);
    edgeX1=reshape(Mi(:,4),4*runInfo.nPtsPerSide-3,[]);

    %------------------------------------------------------------------

    % some of the center points or edge points may extend beyond the image
    % find these and remove them
    centerIdx2remove=find(cX0<1 | cX0>imW | cX1<1 | cX1>imW | cY0<1 | cY0>imL | cY1<1 | cY1>imL);
    [edgeIdx2remove,winIdx2remove]=find(edgeX0<1 | edgeX0>imW | edgeX1<1 | edgeX1>imW | edgeY0<1 | edgeY0>imL | edgeY1<1 | edgeY1>imL);
    idx2remove=unique([centerIdx2remove;winIdx2remove]);

    cX0(idx2remove)=[];
    cX1(idx2remove)=[];
    cY0(idx2remove)=[];
    cY1(idx2remove)=[];

    edgeX0(:,idx2remove)=[];
    edgeX1(:,idx2remove)=[];
    edgeY0(:,idx2remove)=[];
    edgeY1(:,idx2remove)=[];

    nWindows=nWindows-length(idx2remove);

    % set velocity at locations outside cell to zero
    %     ind = sub2ind(size(cellMask0), ceil(cY0), ceil(cX0));
    %     outside = cellMask0(ind) ~= 1;
    %     cY1(outside) = cY0(outside);
    %     cX1(outside) = cX0(outside);
    %
    %     % set velocity at locations outside cell to zero
    %     ind = sub2ind(size(cellMask0), ceil(edgeY0), ceil(edgeX0));
    %     outside = cellMask0(ind) ~= 1;
    %     edgeY1(outside) = edgeY0(outside);
    %     edgeX1(outside) = edgeX0(outside);

    if DEBUG__==1 && pairNumber==1
        % for first frame, show filtered set of vectors in blue
        % show interpolated window centers in red
        % show interpolated window edge points in green
        imshow(vectorCoverageMask); hold on
        quiver(pyxFiltered(:,2),pyxFiltered(:,1),vyxFiltered(:,2),vyxFiltered(:,1),0,'b')
        quiver(cX0,cY0,cX1-cX0,cY1-cY0,0,'r')
        quiver(edgeX0,edgeY0,edgeX1-edgeX0,edgeY1-edgeY0,0,'g')
    end


    %------------------------------------------------------------------

    % Ai: nWin-vector containing area of all the windows in frame i
    % A0=polyarea(edgeX0,edgeY0)';
    % A1=polyarea(edgeX1,edgeY1)';


    % Ii: nWin-vector containing sum of pixel intensities for all windows in
    % frame i
    I0=zeros(nWindows,1);
    I1=zeros(nWindows,1);
    % consider all image pixels for in/out/on test
    [X Y]=meshgrid(1:imW,1:imL);
    for i=1:nWindows
        % get which pixels are inside/on window i's boundary in each frame
        [in0 on0] = inpolygon(X,Y,edgeX0(:,i),edgeY0(:,i));
        [in1 on1] = inpolygon(X,Y,edgeX1(:,i),edgeY1(:,i));

        % we assume that if a pixel center is on a boundary, it is shared by at
        % most 2 adjacent windows. in practice, very few pixel centers will be
        % located on a boundary, and it will be extremely rare that a pixel
        % center forms a vertex for three or more windows. thus, we can sum
        % whole pixel intensities from the "in" list and half pixel intensities
        % (because shared by 2 windows) from the "on" list.

        % first remove "on" (boundary) pixels from "in" list
        in0=logical(in0-on0);
        in1=logical(in1-on1);
        % now sum the intensities
        I0(i)=sum(img0(in0(:)))+sum(img0(on0(:)))/2;
        I1(i)=sum(img1(in1(:)))+sum(img1(on1(:)))/2;
    end

    % voila, it's really that easy
    netTurnover = (I1-I0)./(runInfo.winL.^2);

    %------------------------------------------------------------------

    % get list of indices for windows with net poly or net depoly
    netPolyWinIdx=find(netTurnover>0);
    netDepolyWinIdx=find(netTurnover<=0);

    % for net DEPOLY, we want to record the depoly value of a measurement
    % window where the intensity  last appears, which is the window center in
    % the FIRST frame. for net POLY, we want to record the poly value of a
    % measurement window where the intensity first appears, which is the window
    % center in the SECOND frame.

    % initialize vectors to hold y,x-coordinates of the locations where
    % poly/depoly should be recorded.
    accumY=zeros(size(netTurnover));
    accumX=zeros(size(netTurnover));

    % depoly, so take the first frame centers
    accumY(netDepolyWinIdx)=cY0(netDepolyWinIdx);
    accumX(netDepolyWinIdx)=cX0(netDepolyWinIdx);

    % poly, take location based on flag
    if runInfo.polyLocFlag==0 % take first frame centers
        accumY(netPolyWinIdx)=cY1(netPolyWinIdx);
        accumX(netPolyWinIdx)=cX1(netPolyWinIdx);
    elseif runInfo.polyLocFlag==1 % take the second frame centers
        accumY(netPolyWinIdx)=cY1(netPolyWinIdx);
        accumX(netPolyWinIdx)=cX1(netPolyWinIdx);
    end

    accumTurnover=netTurnover;

    % it might occur that one pixel location has several values associated with
    % it (e.g. poly in window a and depoly in window b map to same pixel). thus
    % we sum these values, since the optical effect of poly/depoly in the same
    % location is additive.  here we find repeated entries, sum their turnover
    % scores into one location, and remove the other locations.
    accumY=round(accumY); accumY(accumY<1)=1; accumY(accumY>imL)=imL;
    accumX=round(accumX); accumX(accumX<1)=1; accumX(accumX>imW)=imW;
    pixIdx=xy2index(accumX,accumY,imL,imW);
    [uniqueE,nOccurences] = countEntries(pixIdx,1,0);
    repeats=uniqueE(nOccurences>1);
    mark4removal=zeros(100000,1); c=1;
    for i=1:length(repeats)
        idxList=find(pixIdx==repeats(i));
        accumTurnover(idxList(1))=sum(accumTurnover(idxList));
        mark4removal(c:c+length(idxList(2:end))-1)=idxList(2:end);
        c=c+length(idxList(2:end));
    end
    mark4removal(mark4removal==0)=[];
    accumY(mark4removal)=[];
    accumX(mark4removal)=[];
    pixIdx(mark4removal)=[];
    accumTurnover(mark4removal)=[];


    
    % interpolate poly/depoly values at every pixel in vectorCoverageMask
    dataM = imDataMap([max(accumY)-min(accumY)+1 max(accumX)-min(accumX)+1],...
        [accumY-min(accumY)+1 accumX-min(accumX)+1],accumTurnover,...
        'infLen',11,...   
        'grid',[5 5],...
        'mask',vectorCoverageMask(min(accumY):max(accumY),min(accumX):max(accumX)));
    % put data back into imL x imW image
    polyDepoly=nan*zeros(imL,imW);
    polyDepoly(min(accumY):max(accumY),min(accumX):max(accumX))=dataM;
   
    runInfo.polyDepolyMovieMin=nanmin(runInfo.polyDepolyMovieMin,nanmin(polyDepoly(:)));
    runInfo.polyDepolyMovieMax=nanmax(runInfo.polyDepolyMovieMax,nanmax(polyDepoly(:)));
    
    if DEBUG__==1
        % plot the data using red/green values from this frame only
        % show outliers in cyan, raw filtered vectors in yellow
        % show cell outlines from frame pair in blue
        cMap=isomorphicColormap('g/r',128);
        cMap=[cMap; [0 0 0]];
        opacity=1;
        [img3C,pixClasses]=imDataMapOverlay(img0,polyDepoly,[-nanmax(abs(polyDepoly(:))),nanmax(abs(polyDepoly(:)))],cMap,opacity);
        img3C(runInfo.edgePix{pairNumber}+2*imL*imW)=1; % show cell from pair i in dark blue
        img3C(runInfo.edgePix{pairNumber+1}+2*imL*imW)=.5; % show cell from pair i+1 in lighter blue
        hmain=get(0,'CurrentFigure');
        if pairNumber==1 && ~isempty(hmain)
            figure(hmain+1);
        else
            figure(1);
        end
        imshow(img3C);
        hold on
        quiver(vecPxyVxy(:,1),vecPxyVxy(:,2),vecPxyVxy(:,3),vecPxyVxy(:,4),0,'c')
        quiver(pyxFiltered(:,2),pyxFiltered(:,1),vyxFiltered(:,2),vyxFiltered(:,1),0,'y')
        % save matlab figure of vector overlay
        saveas(gcf,[runInfo.turnDir filesep 'overlay' filesep 'overlay' indxStr1 '.fig']);
    end

    % save raw data - accumY, accumX, and accumTurnover
    save([runInfo.turnDir filesep 'raw' filesep 'accumY' indxStr1 '.mat'],'accumY');
    save([runInfo.turnDir filesep 'raw' filesep 'accumX' indxStr1 '.mat'],'accumX');    
    save([runInfo.turnDir filesep 'raw' filesep 'accumTurnover' indxStr1 '.mat'],'accumTurnover');
    % save interpolated data
    save([runInfo.turnDir filesep 'interp' filesep 'polyDepoly' indxStr1 '.mat'],'polyDepoly');
    % save vector coverage masks
    imwrite(vectorCoverageMask,[runInfo.turnDir filesep 'vectorCoverageMask' filesep 'vectorCoverageMask' indxStr1 '.tif']);

end
save([runInfo.turnDir filesep 'runInfo'],'runInfo');
save([runInfo.turnDir filesep 'filteredFlow'],'filteredFlow');



% function filteredIm = filterRegion(im, mask, kernel)
% warning('off', 'MATLAB:divideByZero');
% im(~mask) = 0;
% filteredIm = imfilter(im, kernel,'symmetric');
% W = imfilter(double(mask), kernel);
% filteredIm = filteredIm ./ W;
% filteredIm(~mask) = 0;
