function [spy,spx]=spFitVecField(interpVecField)
%SPFITVECFIELD gives spline representation of an interpolated vector field
%
% SYNOPSIS: [spy,spx]=spFitVecField(interpVecField)
%
% INPUT: interpVecField : nVectors x 4 array
%                         [y0 x0 y1 x1] where (y0,x0) is the base and
%                         (y1,x1) is the tip of the vector. A value is
%                         needed at every grid point.
%
% OUTPUT: 
%          spy  : spline fit of y-components of vector field 
%          spx  : spline fit of x-components of vector field
%
% MATLAB VERSION (originally written on): 7.0.1.24704 (R14) Service Pack 1 Windows_NT
%
% USERNAME: kathomps
% DATE: 15-Feb-2006
%
% COMMENTS: currently requires all points on a rectangular grid to be in
% interpVecField

% positions of vector tail
py=interpVecField(:,1);
px=interpVecField(:,2);

% length of vector
vy=interpVecField(:,3)-interpVecField(:,1);
vx=interpVecField(:,4)-interpVecField(:,2);

% get how many grid pts are in y direction
ptsY=length(find(px==px(1))); 

% reshape from n x 2 to grid
gridPY=reshape(py,ptsY,[]);
gridPX=reshape(px,ptsY,[]);
gridVY=reshape(vy,ptsY,[]);
gridVX=reshape(vx,ptsY,[]);

% REPRESENT INTERPOLATED VECTOR FIELD WITH SPLINES
% gridVY(isnan(gridVY)) = 0;
% gridVX(isnan(gridVX)) = 0;
knotsX = optknt(gridPX(1,:),4);
knotsY = optknt(gridPY(:,1)',4);

% least squares spline
spx=spap2({knotsY,knotsX},[4 4],{gridPY(:,1)',gridPX(1,:)},gridVX);
spy=spap2({knotsY,knotsX},[4 4],{gridPY(:,1)',gridPX(1,:)},gridVY);

