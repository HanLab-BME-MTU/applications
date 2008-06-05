function [yTLC,xTLC,yBRC,xBRC,nWinL,nWinW]=findCalcRegion(roiMask,winL)


[imL,imW]=size(roiMask);

% find the limits of the ROI in the mask
[y x]=find(roiMask);
minY=min(y);
maxY=max(y);
minX=min(x);
maxX=max(x);

rangeY=maxY-minY+1;
rangeX=maxX-minX+1;

% dilate the ROI and find the new limits
roiMaskDil=bwmorph(roiMask,'dilate',2*winL);
[y x]=find(roiMaskDil);
minYd=min(y);
maxYd=max(y);
minXd=min(x);
maxXd=max(x);

rangeYd=maxYd-minYd+1;
rangeXd=maxXd-minXd+1;

% how many pixels to add along original ROI to give a large enough
% rectangle around the ROI. the size must contain an integer number of
% windows
maxWinTraveled=3;
halfY1=ceil((maxWinTraveled*winL-mod(rangeY,winL))/2);
halfX1=ceil((maxWinTraveled*winL-mod(rangeX,winL))/2);
halfY2=floor((maxWinTraveled*winL-mod(rangeY,winL))/2);
halfX2=floor((maxWinTraveled*winL-mod(rangeX,winL))/2);

% test
mod(halfY1+halfY2+rangeY,winL);
mod(halfX1+halfX2+rangeX,winL);

% find cushion between circumscribing rectangles around roi and roiDilated
topDiff=minY-minYd;
bottomDiff=maxYd-maxY;
leftDiff=minX-minXd;
rightDiff=maxXd-maxX;

% place ROI in vertical direction
if topDiff>=halfY1 & bottomDiff>=halfY2
    % both fit, no problem
    newMinY=minY-halfY1;
    newMaxY=maxY+halfY2;
elseif topDiff<halfY1 & bottomDiff<halfY2
    % neither fits, just use whole image
    newMinY=floor(mod(imL,winL)/2);
    if newMinY==0
        newMinY=1;
    end
    newMaxY=newMinY+winL*floor(imL/winL)-1;
elseif topDiff<halfY1
    % top doesn't fit, bottom ok
    newMinY=1;
    newMaxY=winL*ceil(maxYd/winL);
    if newMaxY>imL
        newMaxY=newMaxY-winL;
    end
elseif bottomDiff<halfY1
    % bottom doesn't fit, top ok
    newMaxY=imL;
    newMinY=imL-winL*(ceil((imL-minYd)/winL))+1;
    if newMinY<1
        newMinY=newMinY+winL;
    end
end

% place ROI in horizontal direction
if leftDiff>=halfX1 & rightDiff>=halfX2
    % both fit, no problem
    newMinX=minX-halfX1;
    newMaxX=maxX+halfX2;
elseif leftDiff<halfX1 & rightDiff<halfX2
    % neither fits, just use whole image
    newMinX=floor(mod(imW,winL)/2);
    if newMinX==0
        newMinX=1;
    end
    newMaxX=newMinX+winL*floor(imW/winL)-1;
elseif leftDiff<halfX1
    % left doesn't fit, right ok
    newMinX=1;
    newMaxX=winL*ceil(maxXd/winL);
    if newMaxX>imW
        newMaxX=newMaxX-winL;
    end
elseif rightDiff<halfX1
    % right doesn't fit, left ok
    newMaxX=imW;
    newMinX=imW-winL*(ceil((imW-minXd)/winL))+1;
    if newMinX<1
        newMinX=newMinX+winL;
    end
end

% make sure integer number of windows along length and width
nWinL=(newMaxY-newMinY+1)/winL;
nWinW=(newMaxX-newMinX+1)/winL;
if floor(nWinL)~=nWinL || floor(nWinW)~=nWinW
    error('findCalcRegion: non-integer number of windows was calculated')
end
if newMinX<1 || newMinY<1 || newMaxX>imW || newMaxY>imL
    error('findCalcRegion: region extends beyond image boundaries')
end

% yx-coordinates for the top-left and bottom-right corners
yTLC=newMinY; xTLC=newMinX;
yBRC=newMaxY; xBRC=newMaxX;


