function [sqImg,offsetX,offsetY] = makeSquareImage(rectImg,val)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%

% default value for additional image pixel
if nargin < 2
    val=NaN;
end

offsetX=0;
offsetY=0;

% sx: number of rows, sy: number of columns
[sx,sy]=size(rectImg);
sMax=max([sx,sy]);

% remember which side is largest
if sMax == sx
    sxIsMax=true;
else
    sxIsMax=false;
end

sqImg=ones(sMax,sMax)*val;

if sxIsMax
    shift=floor((sMax-sy)/2)+1;
    sqImg(:,shift:shift+sy-1)=rectImg;
    offsetX=shift-1;
else
    shift=floor((sMax-sx)/2)+1;
    sqImg(shift:shift+sx-1,:)=rectImg;
    offsetY=shift-1;
end

if isEven(sMax)
    sqImg(sMax+1,:)=val;
    sqImg(:,sMax+1)=val;
end