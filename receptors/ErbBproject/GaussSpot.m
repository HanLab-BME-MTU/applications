function [ Grid ] = GaussSpot( x, y, s, A)
%GaussSpot creates a 7x7 array with values that represents a symetric 2D gaussian
%with amplitude A sigma s and a sub pixel shift of [x,y]

%Create vectors that cover the corrdinates of x and y
xv = -3:3;
yv = xv;

%shift x and y offsets to center gaussian
x = x-0.5;
y = y-0.5;

%x and y of the guassian are separable.
%Here the integral of the norm distribution is calulated for each pixel in
%x and y
GaussX = normcdf(xv,x,s)-normcdf(xv-1,x,s);
GaussY = normcdf(yv,y,s)-normcdf(yv-1,y,s);

%Finally the two vectors are multiplied into a matirx and scaled by A
Grid = A*(GaussY'*GaussX);

end
