function plotNuclei(nuclei,num)
% Use isosurface to plot cropped nuclei
%   Input:coordinates(x,y,z) and value(v) of each single dot

v=nuclei(num).cube;
[x,y,z]=meshgrid(1:size(v,1),1:size(v,2),1:size(v,3));

figure,
p = patch(isosurface(x,y,z,v));
% why not -3? see doc isosurface
isonormals(x,y,z,v,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud

end

