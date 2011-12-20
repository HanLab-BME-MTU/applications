function [ Grid ] = GaussSpotGen( x, y, xs, ys, A, pix)
%GaussSpotGen creates a pix x pix array with values that represents an
%asymetric 2D gaussian with amplitude A sigma s and a sub pixel shift of
%[x,y]. Pix must be a odd number

%defines the range for grid size
pN = floor(pix/2);

%Create vectors that cover the corrdinates of x and y
xv = -pN:pN;
yv = xv;

%shift x and y offsets to center gaussian
x = x-0.5;
y = y-0.5;

%x and y of the guassian are separable.
%Here the integral of the norm distribution is calulated for each pixel in
%x and y
GaussX = normcdf(xv,x,xs)-normcdf(xv-1,x,xs);
GaussY = normcdf(yv,y,ys)-normcdf(yv-1,y,ys);

%Finally the two vectors are multiplied into a matirx and scaled by A
Grid = A*(GaussY'*GaussX);

end
