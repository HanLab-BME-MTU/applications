function [pos_f,force]=reg_FEM2D_TFM(E,v,grid_mat,iu_mat,i_max,j_max)
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
structModel2D=createpde('structural','static-solid');

% //Define geometry *******************************************************
length = 198;
width = 198;
R = [3 4 -length length length -length ...
          -width -width  width   width]';
g = decsg(R,'R','R');
geometryFromEdges(structModel2D,g);

% //Specify material properties *******************************************
structuralProperties(structModel2D,'YoungsModulus',E,'PoissonsRatio',v);

% //Apply fixed BC on bottom face *****************************************
structuralBC(structModel2D,'Edge',3,'XDisplacement',0);

% //Apply displacement field on top face **********************************

%Change xLoc and yLoc to come from grid_mat
[xLoc,yLoc]=ndgrid(-198/2:198/2,-198/2:198/2);
%Add disp_x and disp_y as extracted from iu_mat
disp_x = ux(:,:,1);
disp_y = uy(:,:,1);
dispInterpX=griddedInterpolant(xLoc,yLoc,disp_x);
dispInterpY=griddedInterpolant(xLoc,yLoc,disp_y);

%define function handle
hydrogelDisp = @(location,state)[dispInterpX(location.x,location.y); ...
                                 dispInterpY(location.x,location.y);];
                                 
%pass function handle and define BCs
structuralBoundaryLoad(structModel2D,'Face',1,'Displacement',hydrogelDisp,'Vectorize','on');

% //Generate mesh for model ***********************************************
generateMesh(structModel2D,'Hmax',10);

% //Solving structural model **********************************************
structModel2DResults=solve(structModel2D);

% //Visualizing results ***************************************************
figure(2)
pdeplot(structModel2D,'ColorMapData',structModel2DResults.Displacement.Magnitude, ...
    'Deformation',structModel2DResults.Displacement)
end