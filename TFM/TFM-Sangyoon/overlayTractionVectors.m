function [] = overlayTractionVectors( forceField, displField,vectorScale,vectorSparsity,varargin)
%function [] = overlayTractionVectors(forceField, vectorScale,
% vectorSparsity) overlays vectors preferrably on top of traction map. 
% 
% input:    forceField:                 forceField with fields containing
%                                       pos and vec
%           displField:                 displacement field that will be
%                                       used to calculate the true position
%                                       of force sources. 
%           vectorScale:                vector display scale (default=1)
%           vectorSparsity:             Gap between the vectors (default=2)
%           varargin:                   vector color (defalut=purple)
% output:   
%           quiver plot on existing, or newly generated figure
% 
% It assumes traction map is already displayed. So it shifts the vector
% field based on the upper left position of the grid in the forceField.
% 
% Sangyoon Han, March 2017

% Input
if nargin<4
    vectorSparsity=2;
end
if nargin<3
    vectorSparsity=2;
    vectorScale = 1;
end
ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Color',[75/255 0/255 130/255]);
ip.addParamValue('ShiftField',true,@islogical);
ip.parse(varargin{:});
vecColor=ip.Results.Color;
shiftField=ip.Results.ShiftField;

% if isempty(varargin)
%     vecColor=[75/255 0/255 130/255];
% else
%     vecColor=varargin{2};
% end
reg_grid =createRegGridFromDisplField(displField,4); %2=2 times fine interpolation

% Make the detailed map of force and displacement
[grid_mat,~,~,~] = interp_vec2grid(displField.pos, displField.vec,[],reg_grid);    
    
reg_grid_coarse =createRegGridFromDisplField(displField,1/vectorSparsity); %2=2 times fine interpolation
[grid_mat_coarse,iu_mat_coarse,~,~] = interp_vec2grid(displField.pos, displField.vec,[],reg_grid_coarse);
pos_coarse = [reshape(grid_mat_coarse(:,:,1),[],1) reshape(grid_mat_coarse(:,:,2),[],1)]; %dense
disp_vec_coarse = [reshape(iu_mat_coarse(:,:,1),[],1) reshape(iu_mat_coarse(:,:,2),[],1)]; 
[~,if_mat_coarse,~,~] = interp_vec2grid(forceField.pos, forceField.vec,[],reg_grid_coarse);
force_vec_coarse = [reshape(if_mat_coarse(:,:,1),[],1) reshape(if_mat_coarse(:,:,2),[],1)]; 

[~,tmat_coarse, ~, ~] = interp_vec2grid(pos_coarse+disp_vec_coarse, force_vec_coarse,[],grid_mat_coarse); %1:cluster size

tmat_vecx = reshape(tmat_coarse(:,:,1),[],1);
tmat_vecy = reshape(tmat_coarse(:,:,2),[],1);
pos_vecx = reshape(grid_mat_coarse(:,:,1),[],1);
pos_vecy = reshape(grid_mat_coarse(:,:,2),[],1);
% forceScale=0.1*max(sqrt(tmat_vecx.^2+tmat_vecy.^2));
% quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), vectorScale*tmat_vecx./forceScale,...
%     vectorScale*tmat_vecy./forceScale,0,'Color',vecColor);
if ~shiftField
    grid_mat(1,1,1)=0;
    grid_mat(1,1,2)=0;
end
quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), vectorScale*tmat_vecx,...
    vectorScale*tmat_vecy,0,'Color',vecColor);
%hq = % hq.ShowArrowHead = 'off';
    
return;
