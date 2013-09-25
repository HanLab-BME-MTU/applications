function [A,APer] = meshSurfaceArea(f,v)
%MESHSURFACEAREA calculates the surface area of the input triangulated mesh
%
% area = meshSurfaceArea(faces,vertices)
% [area,areaPerFace] = meshSurfaceArea(faces,vertices)
%
% This just calculates the area of each triangle on the mesh and sums them
% up. (e.g. the type of mesh that is output by isosurface.m etc.)
%
%   Input:
%
%       vertices - the XYZ coordinates of the vertices of the mesh. 
%   
%       faces - indices specifying the vertices which make up each
%       triangular face of the mesh.
%
%   Output:
%
%       area - Total area of all faces (surface area)
%
%       areaPerFace - area of each face.
%
% Hunter Elliott
% 4/2013


%Get vector AB connecting first two vertices
ab = v(f(:,2),:) - v(f(:,1),:);
%..and AC connecting first and third
ac = v(f(:,3),:) - v(f(:,1),:);

%Magnitude of cross is area of parallelogram, area of triangle is 1/2 that:
APer = 1/2 * sqrt(sum(cross(ab,ac,2) .^2 ,2));
A = sum(APer);


