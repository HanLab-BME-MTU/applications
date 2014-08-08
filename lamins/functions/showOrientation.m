function [ thetaDegrees, cm ] = showSteerableFilterOrientation( theta, nms, mask )
%showOrientation Show the orientation of the steerable filter projected
%onto the nms

if(nargin < 3)
    mask = ones(size(theta));
end

% convert to degrees. Should be in the range (-90,90]
thetaDegress = theta*180/pi;
% set masked portions to -90 degrees
thetaDegrees(nms > 0 & mask) = -90;

% show the image with clim -90 to 90
imshow(thetaDegrees,[-90 90]);
% 
cm = [ [0 0 0]; hsv(180)];
colormap(cm);


end

