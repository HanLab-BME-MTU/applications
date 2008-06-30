function [polyDepoly,accumY,accumX,runInfo]=turnoverMap(runInfo,nImgs2Analyze,winL,methodStr)


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

% check input of methodStr
if nargin<4 || isempty(methodStr)
    runInfo.methodStr='frame0grid'; % default - will plot on original grid
else
    runInfo.methodStr=methodStr;
end

% normalized image directory
normDir=[runInfo.imDir filesep 'norm' filesep 'normMats'];

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
mkdir([runInfo.turnDir filesep 'raw']);
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

%------------------------------------------------------------------
% ITERATE THROUGH FRAME PAIRS
for pairNumber=1:runInfo.nImgs2Analyze-1
    
    indxStr1=sprintf(strg,pairNumber);   
    
    %------------------------------------------------------------------
    % READ CELL MASKS AND LOAD NORMALIZED IMAGES FOR FRAMES i AND i+1

    if pairNumber==1
        % mask0: dilated polyDepolyROI x vectorCoverageMask
        fileNameMask0=[char(listOfVecCovMasks(pairNumber,2)) filesep char(listOfVecCovMasks(pairNumber,1))];
        mask0=double(imread(fileNameMask0)).*roiMaskDil;
        
        % load first frame image
        indxStr1=sprintf(strg,pairNumber);
        img0=load([normDir filesep 'norm_image' indxStr1]);
        img0=img0.normImg;
    else
        mask0=mask1;
        img0=img1;
    end
    % mask1: dilated polyDepolyROI x vectorCoverageMask
    fileNameMask1=[char(listOfVecCovMasks(pairNumber+1,2)) filesep char(listOfVecCovMasks(pairNumber+1,1))];
    mask1=double(imread(fileNameMask1)).*roiMaskDil;
    
    % load second frame image
    indxStr2=sprintf(strg,pairNumber+1);
    img1=load([normDir filesep 'norm_image' indxStr2]);
    img1=img1.normImg;

    % smooth the images
    blurKernel=fspecial('gaussian', 11, 5);
    img0 = imfilter(img0, blurKernel,'symmetric');
    img1 = imfilter(img1, blurKernel,'symmetric');

    %------------------------------------------------------------------


    
    % get grid for calculating poly/depoly
    % findCalcRegion gives coordinates for a rectangle surrounding the
    % vectorCoverageMasks, the sides of which are a multiple of gridSpace.
    % we need to add 1 to each side to allow for integer coordinates. for
    % example, if the rectangle were 15 pixels long, and gridSpace=5; we
    % want to end up with coordinates at 1,6,11,16.
    gridSpace=5;
    [yTLC,xTLC,yBRC,xBRC,nWinL,nWinW]=findCalcRegion(mask0 | mask1,gridSpace);
    if yBRC<runInfo.imL
        yBRC=yBRC+1;
    elseif yTLC>1
        yTLC=yTLC-1;
    else
        warning('TURNOVERMAP: border problem in grid to calculate poly/depoly');
    end

    if xBRC<runInfo.imW
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
    pyxFiltered=filteredFlow(1,pairNumber).pyxFilt;
    vyxFiltered=filteredFlow(1,pairNumber).vyxFilt;    
    
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
    
    % how many points on each side of the measurement window to use
    runInfo.nPtsPerSide=4; 
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
    centerIdx2remove=find(cX0<1 | cX0>runInfo.imW | cX1<1 | cX1>runInfo.imW | cY0<1 | cY0>runInfo.imL | cY1<1 | cY1>runInfo.imL);
    [edgeIdx2remove,winIdx2remove]=find(edgeX0<1 | edgeX0>runInfo.imW | edgeX1<1 | edgeX1>runInfo.imW | edgeY0<1 | edgeY0>runInfo.imL | edgeY1<1 | edgeY1>runInfo.imL);
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
    [X Y]=meshgrid(1:runInfo.imW,1:runInfo.imL);
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
    if isequal(runInfo.methodStr,'frame0grid') % take first frame centers
        accumY(netPolyWinIdx)=cY0(netPolyWinIdx);
        accumX(netPolyWinIdx)=cX0(netPolyWinIdx);
    elseif isequal(runInfo.methodStr,'frame1grid') % take the second frame centers
        accumY(netPolyWinIdx)=cY1(netPolyWinIdx);
        accumX(netPolyWinIdx)=cX1(netPolyWinIdx);
    end

    accumTurnover=netTurnover;

    % it might occur that one pixel location has several values associated with
    % it (e.g. poly in window a and depoly in window b map to same pixel). thus
    % we sum these values, since the optical effect of poly/depoly in the same
    % location is additive.  here we find repeated entries, sum their turnover
    % scores into one location, and remove the other locations.
    accumY=round(accumY); accumY(accumY<1)=1; accumY(accumY>runInfo.imL)=runInfo.imL;
    accumX=round(accumX); accumX(accumX<1)=1; accumX(accumX>runInfo.imW)=runInfo.imW;
    pixIdx=xy2index(accumX,accumY,runInfo.imL,runInfo.imW);
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


    maskBoth=mask0 | mask1;
    % interpolate poly/depoly values at every pixel in vectorCoverageMask
    dataM = imDataMap([max(accumY)-min(accumY)+1 max(accumX)-min(accumX)+1],...
        [accumY-min(accumY)+1 accumX-min(accumX)+1],accumTurnover,...
        'infLen',11,...   
        'grid',[5 5],...
        'mask',maskBoth(min(accumY):max(accumY),min(accumX):max(accumX)));
    % put data back into imL x imW image
    polyDepoly=nan*zeros(runInfo.imL,runInfo.imW);
    polyDepoly(min(accumY):max(accumY),min(accumX):max(accumX))=dataM;
   
    runInfo.polyDepolyMovieMin=nanmin(runInfo.polyDepolyMovieMin,nanmin(polyDepoly(:)));
    runInfo.polyDepolyMovieMax=nanmax(runInfo.polyDepolyMovieMax,nanmax(polyDepoly(:)));
  

    % save raw data - accumY, accumX, and accumTurnover
    save([runInfo.turnDir filesep 'raw' filesep 'accumY' indxStr1 '.mat'],'accumY');
    save([runInfo.turnDir filesep 'raw' filesep 'accumX' indxStr1 '.mat'],'accumX');    
    save([runInfo.turnDir filesep 'raw' filesep 'accumTurnover' indxStr1 '.mat'],'accumTurnover');
    % save interpolated data
    save([runInfo.turnDir filesep 'interp' filesep 'polyDepoly' indxStr1 '.mat'],'polyDepoly');

end
save([runInfo.turnDir filesep 'runInfo'],'runInfo');

