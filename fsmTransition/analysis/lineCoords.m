function P=lineCoords(A,B,sampling)
% lineCoords returns all pixel coordinates along a line defined by two points A and B
%
% P=lineCoords(A,B,sampling)
%
% INPUT
%
% A         : origin point coordinates A=[y x ...], any number of dimensions
% A         : origin point coordinates B=[y x ...], anynumber of dimensions
% sampling  : (pixels) distance between two adjacent points along the profile
%             default: 1
%
% OUTPUT
% P         : coordinates of the points along the profile
%                         [y1 y2 y3 ... yn;
%                     P =  x1 x2 x3 ... xn;
%                                ...      ] 
% Aaron Ponti, 12/08/2003

% Check input
if nargin==2
    sampling=1;
end

if size(A,2)>size(A,1)
    A=A';
end

if size(B,2)>size(B,1)
    B=B';
end

if size(A)~=size(B)
    error('The points must have the same number of coordinates.');
end

% Numer of dimensions
dim=size(A,1);

% Calculate vector
V=B-A;

% Calculate vector length
lenV=sqrt(dot(V,V));

% Sample
steps=([0:sampling/lenV:1])';
len=length(steps);

% Calculate pixel coordinates along line
for i=1:dim
    P(1:len,i)=A(i)+steps.*V(i);
end
