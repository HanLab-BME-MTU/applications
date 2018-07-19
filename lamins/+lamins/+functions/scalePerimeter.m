function [ scaled_bw, scaledpos ] = scalePerimeter( bw, scale )
%scalePerimeter Scales the perimeter of the object given in bw

boundary = bwboundaries(bw);
boundary = boundary{1};
centroid = mean(boundary);

% translate
centroid = repmat(centroid,length(boundary),1);
boundary = boundary - centroid;
[theta,rho] = cart2pol(boundary(:,1),boundary(:,2));
[scaledpos(:,1), scaledpos(:,2)] = pol2cart(theta,rho*scale);
scaledpos = scaledpos + centroid;
scaled_bw = zeros(size(bw));
% filter = scaledpos(:,1) < 
scaled_bw(sub2ind(size(bw),round(scaledpos(:,1)),round(scaledpos(:,2)))) = 1;


end

