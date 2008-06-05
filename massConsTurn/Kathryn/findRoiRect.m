function [yTLC,xTLC,yBRC,xBRC]=findRoiRect(roiMask,winL)
%FINDROIRECT returns corner coordinates of ROI-circumscribing rectangle
%
% DESCRIPTION: findRoiRect finds the best rectangle to circumscribe a
% region of interest (ROI) in a binary mask. The rectangle's sides must be
% divisible by winL.
%
% SYNOPSIS: [yTLC,xTLC,yBRC,xBRC]=findRoiRect(roiMask,winL)
%
% INPUT: roiMask   : binary mask containing one blob of 1's which
%                    represent something like a cell
%        winL      : integer by which the rectangle's dimensions must
%                    divide evenly
%
% OUTPUT: yTLC,xTLC: yx-coordinates of Top Left Corner of the rectangle
%         yBRC,xBRC: yx-coordinates of Bottom Right Corner of the rectangle
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
% USERNAME: kathomps
% DATE: 13-Jun-2007

roiMask=mask0;
[imL,imW]=size(roiMask);

% find the limits of the ROI in the mask
[y x]=find(roiMask);
minY=min(y);
maxY=max(y);
minX=min(x);
maxX=max(x);

rangeY=maxY-minY+1;
rangeX=maxX-minX+1;

% how many pixels to add in x and y to make integer number of boxes with a
% wide enough border around the ROI
amt2addY=2*winL-mod(rangeY,winL);
amt2addX=2*winL-mod(rangeX,winL);

% these are the number of pixels to work with in all four directions about
% the ROI
nPixAbove=minY-1;
nPixBelow=imL-maxY;
nPixLeft=minX-1;
nPixRight=imW-maxX;

halfY=ceil(amt2addY/2);
halfX=ceil(amt2addX/2);
topFlag=0; bottomFlag=0; leftFlag=0; rightFlag=0;
if nPixAbove<=halfY
    topFlag=1;
    newMinY=1;
end
if nPixBelow<=halfY
    bottomFlag=1;
    newMaxY=imL;
end
if nPixLeft<=halfX
    leftFlag=1;
    newMinX=1;
end
if nPixRight<=halfX
    rightFlag=1;
    newMaxY=imW;
end

if topFlag==1
    % we set the new minY at 1, so calc new maxY
    newMaxY=


    if imL-rangeY>=amt2addY
        halfY=ceil(amt2addY/2);
        if nPixAbove>=halfY && nPixBelow>=(amt2addY-halfY) % both ok to add half
            minY=minY-halfY;
            maxY=maxY+(amt2addY-halfY);
        elseif nPixAbove<halfY % nPixAbove too small to add half
            minY=1;
            maxY=rangeY+amt2addY;
        else % nPixBelow too small to add half
            maxY=imL;
            minY=imL-(rangeY+amt2addY)+1;
        end
    else % if the amount to add is greater than the available pixels
        minY=1;
        maxY=imL;
    end

    if imW-rangeX>=amt2addX
        halfX=ceil(amt2addX/2);
        if nPixLeft>=halfX && nPixRight>=(amt2addX-halfX) % both ok to add half
            minX=minX-halfX;
            maxX=maxX+(amt2addX-halfX);
        elseif nPixLeft<halfX % nPixLeft too small to add half
            minX=1;
            maxX=rangeX+amt2addX;
        else % nPixRight too small to add half
            maxX=imW;
            minX=imW-(rangeX+amt2addX)+1;
        end
    else % if the amount to add is greater than the available pixels
        minX=1;
        maxX=imW;
    end


    % if it has to be the whole length, cut down to number of boxes that fit
    if minY==1 && maxY==imL
        newRangeY=winL*floor(imL/winL);
        minY=floor((imL-newRangeY)/2);
        maxY=minY+newRangeY-1;
    end

    % if it has to be the whole width, cut down to number of boxes that fit
    if minX==1 && maxX==imW
        newRangeX=winL*floor(imW/winL);
        minX=floor((imW-newRangeX)/2);
        maxX=minX+newRangeX-1;
    end

    % yx-coordinates for the top-left and bottom-right corners
    yTLC=minY; xTLC=minX;
    yBRC=maxY; xBRC=maxX;

