function [ Map ] = GaussField( n, s, A , STN)
%GaussField generates a random discrete "image" with randomly placed
%gaussian functions with guassianly distributed noise
%   Creates an artifical microscope image of n diffraction limited points.
%   The points are modeled as 2D Gaussians with sigmas = s (in pixels).
%   Each points is allowed subpixel positional accuracy and to create the
%   image each pixel contains the integrated intensities of all gaussians
%   that are with in 5*s pixels of the point's center. Each point
%   represented as = A*gaussian((x-center)/s)
%
%   STN,signal to noise determines the sigma of the rayleigh distributed
%   noise

dbstop if all error

%Noise sigma
nS = A/STN;

% Output with noise already added
Map = abs(normrnd(0, nS, 256,256));

%All spots will be completely within the image
List = unifrnd(3.5,252.5,n,2); 

for p = 1:1:n 
    x = List(p,:);
    temp = GaussSpot(mod(x(1),1),mod(x(2),1),s,A);
    
    Map((-3:3)+floor(x(1)),(-3:3)+floor(x(2))) = ...
        Map((-3:3)+floor(x(1)),(-3:3)+floor(x(2)))+temp;
    if mod(p,5) == 0
        b=1;
    end
    imshow(Map, [min(min(Map)) max(max(Map))] );
end

