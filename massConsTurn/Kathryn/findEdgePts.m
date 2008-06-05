function [py,px]=findEdgePts(y0,x0,winL,nPtsPerSide)
% for square centers (y0,x0), get coordinates on square perimeter
%
% INPUT:
%   y0(x0)     : n-vector containing y(x)-coordinates of the center points
%                of n squares 
%   winL       : n-vector containing length of each square's side
%                (or a constant, if all the squares have a uniform size)
%   nPtsPerSide: the number of points to generate per side of the square
%                total number of uniqe points will be 4(nPtsPerSide-1)
%
% OUTPUT:
%   py(px)    : 4(nPtsPerSide-1) x nSquares array containing
%               y(x)-coordinates of the polygon that surrounds the square


% top left corner
tlY=y0-winL./2;
tlX=x0-winL./2;
% top right corner
trY=y0-winL./2;
trX=x0+winL./2;
% bottom right corner
brY=y0+winL./2;
brX=x0+winL./2;
% bottom left corner
blY=y0+winL./2;
blX=x0-winL./2;
% initialize the list of y-coordinates along the four sides for each square
sly=zeros(nPtsPerSide,length(y0));
s2y=zeros(nPtsPerSide,length(y0));
s3y=zeros(nPtsPerSide,length(y0));
s4y=zeros(nPtsPerSide,length(y0));
% iterate through the squares, generating points along the four sides
% initialize the list of x-coordinates along the four sides for each square
slx=zeros(nPtsPerSide,length(y0));
s2x=zeros(nPtsPerSide,length(y0));
s3x=zeros(nPtsPerSide,length(y0));
s4x=zeros(nPtsPerSide,length(y0));
for i=1:length(y0) % iterate through the squares
    s1y(:,i)=linspace(tlY(i),trY(i),nPtsPerSide)'; % top
    s2y(:,i)=linspace(trY(i),brY(i),nPtsPerSide)'; % right
    s3y(:,i)=linspace(brY(i),blY(i),nPtsPerSide)'; % bottom
    s4y(:,i)=linspace(blY(i),tlY(i),nPtsPerSide)'; % left
    
    s1x(:,i)=linspace(tlX(i),trX(i),nPtsPerSide)'; % top
    s2x(:,i)=linspace(trX(i),brX(i),nPtsPerSide)'; % right
    s3x(:,i)=linspace(brX(i),blX(i),nPtsPerSide)'; % bottom
    s4x(:,i)=linspace(blX(i),tlX(i),nPtsPerSide)'; % left
end
% remove repeats at the corners
s2x(1,:)=[]; s2y(1,:)=[];
s3x(1,:)=[]; s3y(1,:)=[];
s4x(1,:)=[]; s4y(1,:)=[];
% concatenate y and x coordinates
py=[s1y; s2y; s3y; s4y];
px=[s1x; s2x; s3x; s4x];



