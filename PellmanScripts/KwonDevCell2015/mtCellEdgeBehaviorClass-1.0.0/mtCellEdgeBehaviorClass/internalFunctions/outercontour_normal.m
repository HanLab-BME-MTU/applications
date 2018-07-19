function [contourYX, normalYX] = outercontour_normal(bwIn)

% in case the input binary image has holes
bwIn = imfill(bwIn, 'holes');

% get the outer boundary of the region
RoiYX = bwboundaries(bwIn);
centerYX = regionprops(bwIn,'Centroid');
contourYX = RoiYX{1};

contour_length = size(contourYX,1);

%Initialize the normal 
normalYX = contourYX*0;

% a parameter to tune, for how much is the neiboring radius (in points) to be considered for direction
r = 10;
contourYX_replicate = [contourYX(end-r+1:end,:); contourYX; contourYX(1:r,:)];

for i = r + 1 : contour_length + r
    neibor_ps = contourYX_replicate(i-r:i+r, :);
    [pc,score] = princomp(neibor_ps);
    normalYX(i-r,:) = [pc(1,2) pc(2,2)];
end % for every point on the contour, find the normal

% for every point on the contour, correct the normal pointing outward using the center to point direction
for i = 1 : contour_length
    if(dot(normalYX(i,:),contourYX(i,:)-[ centerYX.Centroid(2) centerYX.Centroid(1)])<0)
        normalYX(i,:) = -normalYX(i,:);
    end
end 

% for every point on the contour, correct the singularity of the normal
% using the neiboring normals
normal_angle=[];
for i = 1 : contour_length
    normal_angle(i) = atan2(normalYX(i,1),normalYX(i,2));
    if (i>1 && i <contour_length)
        if abs(normal_angle(i) - normal_angle(i-1)) >2.8 && abs(normal_angle(i) - normal_angle(i-1)) < 6
            normalYX(i,:) = -normalYX(i,:);
            normal_angle(i) = atan2(normalYX(i,1),normalYX(i,2));
        end
    end
end 

% But they could all pointing inward, so check for the mean orientation
direct_sign=[];
for i = 1 : contour_length
    direct_sign(i) = sign(dot(normalYX(i,:),contourYX(i,:)-[ centerYX.Centroid(2) centerYX.Centroid(1)]));    
end 

if mean(double(direct_sign))<0
    normalYX = -normalYX;
end







