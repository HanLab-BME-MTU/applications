%Creates a Gaus style heat map from a point process
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 9/6/2011
%
%
%

function [map]=PointP_GausMap(pos,thres)
%PointP_GausMap generates a Gaus Lab style heat map of the point process P
%this means calculating the l(r)-r where r is the thres (distance) which is
%interpolated into a surfacemap.

lr = PointP_Lr(pos,thres)-thres;

xMin = min(pos(:,1));
xMax= max(pos(:,1));
yMin=min(pos(:,2));
yMax = max(pos(:,2));

[XI,YI]= meshgrid(xMin:.2:xMax,yMin:.2:yMax);
map = griddata(pos(:,1),pos(:,2),lr,XI,YI,'v4');

end
