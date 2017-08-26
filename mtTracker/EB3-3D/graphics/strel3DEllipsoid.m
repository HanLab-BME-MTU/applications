function [se,m]=strel3DEllipsoid(lat,ax)
sw = (lat-1)/2;
sn = (ax-1)/2;
ses2 = ceil(lat/2);
sen2 = ceil(ax/2);
[y,x,z] = meshgrid(-sw:sw, -sw:sw, -sn:sn);
m = sqrt((x/sw).^2 + (y/sw).^2 + (z/sn).^2);
b = (m <= m(ses2, lat, sen2));
se = strel('arbitrary', b);