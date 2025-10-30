function []=test_FEM_TFM(ux,uy,E,v)
%reg_FEM_TFM performs FEA using material properties and substrate geometry 
%coupled with measured displacement field to produce estimated force field
%for the given substrate deformation.
%
%   Inputs:
%       E = Young's Modulus of substrate
%       v = Poisson's Ratio of substrate
%       grid_mat = grid of points associated with the displacements in iu_mat
%       iu_mat = matrix of displacement magnitudes
%       i_max = 
%       j_max = 
%
%   Outputs:
%       pos_f = nodal locations
%       force = estimated force at corresonding nodal locations

% //Generate model container **********************************************
structModel=createpde('structural','static-solid');

% //Define geometry *******************************************************
%reduce to micron size geometry, 72 nm/pix
%w = length(grid_mat(:,1));
w = 198;
gm=multicuboid(w,w,15);
structModel.Geometry=gm;

% //Specify material properties *******************************************
structuralProperties(structModel,'YoungsModulus',E,'PoissonsRatio',v);

% //Apply fixed BC on bottom face *****************************************
structuralBC(structModel,'Face',1,'Constraint','fixed');

% //Apply displacement field on top face **********************************

%Change xLoc and yLoc to come from grid_mat
[xLoc,yLoc]=ndgrid(-w/2:w/2-1,-w/2:w/2-1);
%Add disp_x and disp_y as extracted from iu_mat

dispInterpX=griddedInterpolant(xLoc,yLoc,ux(:,:,1));
dispInterpY=griddedInterpolant(xLoc,yLoc,uy(:,:,1));
dispInterpZ=griddedInterpolant(xLoc,yLoc,zeros(198,198));

%define function handle
hydrogelDisp = @(location,state)[dispInterpX(location.x,location.y); ...
                                 dispInterpY(location.x,location.y); ...
                                 dispInterpZ(location.x,location.y)];
                                 
%pass function handle and define BCs
structuralBC(structModel,'Face',2,'Displacement',hydrogelDisp,'Vectorize','on');

% //Generate mesh for model ***********************************************
generateMesh(structModel,'Hmax',10);

% //Solving structural model **********************************************
structModelResults=solve(structModel);

% //Visualizing results ***************************************************
figure(2)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.Magnitude, ...
    'Deformation',structModelResults.Displacement)
end