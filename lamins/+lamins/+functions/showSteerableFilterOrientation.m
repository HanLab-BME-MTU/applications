function [ thetaDegrees, cm ] = showSteerableFilterOrientation( theta, nms, mask, arrows )
%showOrientation Show the orientation of the steerable filter projected
%onto the nms with an optional mask
%
% @param theta is the orientation of the steerable filter and comes from [res, theta, nms ] = steerableDetector
% @param nms is the nonmaximum suppression and comes from [res, theta, nms ] = steerableDetector
% @param mask is a binary logical mask for which false pixels are excluded,
%             use [] to include everything
% @param arrows 0 to turn off all arrows (default), 1 includes all arrows.
%        higher values for arrows reduce the number of arrows plotted
%
% @return thetaDegrees theta in degrees with nms and mask applied, -90
% indicates no value
%
% @cm colormap used to display thetaDegrees, black prepended to hsv(180)

if(nargin < 3 || isempty(mask))
    % accept all pixels
    mask = ones(size(theta));
end
if(nargin < 4)
    % turn off arrows by default
    arrows = 0;
end

% convert to degrees. Should be in the range (-90,90]
thetaDegrees = theta*180/pi;
% set masked portions to -90 degrees
thetaDegrees(nms == 0 | ~mask) = -90;

% show the image with clim -90 to 90
imshow(thetaDegrees,[-90 90]);
% set -90 to black
cm = [ [0 0 0]; hsv(180)];
colormap(cm);

% setup and label colorbar
hcb = colorbar;
set(hcb,'YTick',-90:15:90)
set(get(hcb,'YLabel'),'String','Degrees');

% add arrows
if(arrows > 0)
    % select all pixels remaining after nms and within the mask
    selected = nms ~= 0 & mask;
    idx =  find(selected);
    % grab a random permutation of arrows, reduce the value of arrows
    idx = idx(randperm(length(idx),floor(length(idx)/arrows)));
    % grab the subscripts for the coordinates
    [r,c] = ind2sub(size(theta),idx);
    hold on;
    % xy are reversed in ij space
    % theta is the normal vector
    quiver(c,r,sin(theta(idx)),-cos(theta(idx)),0,'color','white');
end

end

