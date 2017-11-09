function scrable_output_feature = ...
    scrable_network_analysis(VIF_current_model,VIF_current_seg,CellROI, radius,T_sigma,O_sigma)
% function for density feature for scrabled network
% Liya Ding 06.2014.

img_size = size(VIF_current_seg);

if(isempty(CellROI))
    CellROI=ones(img_size);
end

[scrable_digital_model,scrable_orientation_model,scrable_VIF_current_seg,scrable_VIF_current_orientation,...
    scrable_XX,scrable_YY,scrable_OO, scrable_II] ...
    = filament_model_scrable(VIF_current_model,img_size, T_sigma,O_sigma,CellROI);

% scrable_VIF_current_seg = scrable_VIF_current_seg.*CellROI;

density_H = double((fspecial('disk',radius))>0);
density_H = density_H./(sum(sum(density_H)));

density_filament = imfilter(scrable_VIF_current_seg,density_H, 'replicate','same');

scrable_output_feature=[];

scrable_output_feature.density_filament=density_filament;
scrable_output_feature.current_seg=scrable_VIF_current_seg;

