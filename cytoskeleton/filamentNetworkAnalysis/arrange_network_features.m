function output_feature_reorganized = arrange_network_features(output_feature)

%  straightness_per_filament_pool: [1x1611 double]
%                  curvature_per_pixel_pool: [43174x1 double]
%                  length_per_filament_pool: [1x1611 double]
%            pixel_number_per_filament_pool: [1x1611 double]
%                          density_filament: [2160x2560 double]
%                 scrabled_density_filament: [2160x2560 double]
%            orientation_pixel_pool_display: [43174x1 double]
%     orientation_pixel_pool_display_center: [43174x1 double]
%                   filament_mean_curvature: [1x1611 double]
%               intensity_per_filament_pool: [1x1611 double]
%          mean_intensity_per_filament_pool: [1x1611 double]
%           intensity_per_fat_filament_pool: [1x1611 double]
%      mean_intensity_per_fat_filament_pool: [1x1611 double]
%                   scale_per_filament_pool: [1x1611 double]
%                      st_per_filament_pool: [1x1611 double]
%                 mean_st_per_filament_pool: [1x1611 double]
%                  st_per_fat_filament_pool: [1x1611 double]
%             mean_st_per_fat_filament_pool: [1x1611 double]
%                         number_of_nucleus: 30
%                               profileCell: {1x30 cell}
%                            profileAllCell: [1x1 struct]
%                          Centripetal_fila: [1x1611 double]
%                         Centripetal_pixel: [41563x1 double]
%                                 Cell_Mask: [2160x2560 double]



% if there is no curvature calculated, add empty field
if (~isfield(output_feature, 'curvature_per_pixel_pool'))
    output_feature.curvature_per_pixel_pool  = [];
end

% if there is no profileCell calculated, add empty field
if (~isfield(output_feature, 'profileCell'))
    output_feature.profileCell=[];
end

% if there is no Centripetal_fila calculated, add empty field
if (~isfield(output_feature, 'Centripetal_fila'))
    output_feature.Centripetal_fila  = [];
end

% if there is no number of nucleus calculated, add empty field
if (~isfield(output_feature, 'number_of_nucleus'))
    output_feature.number_of_nucleus  = [];
end

if (~isfield(output_feature, 'straightness_per_filament_pool'))
    output_feature.straightness_per_filament_pool  = [];
end


if (~isfield(output_feature, 'length_per_filament_pool'))
    output_feature.length_per_filament_pool  = [];
end

if (~isfield(output_feature, 'pixel_number_per_filament_pool'))
    output_feature.pixel_number_per_filament_pool  = [];
end

if (~isfield(output_feature, 'intensity_per_filament_pool'))
    output_feature.intensity_per_filament_pool  = [];
end

if (~isfield(output_feature, 'intensity_per_fat_filament_pool'))
    output_feature.intensity_per_fat_filament_pool  = [];
end

if (~isfield(output_feature, 'scale_per_filament_pool'))
    output_feature.scale_per_filament_pool  = [];
end

if (~isfield(output_feature, 'st_per_filament_pool'))
    output_feature.st_per_filament_pool  = [];
end

if (~isfield(output_feature, 'st_per_fat_filament_pool'))
    output_feature.st_per_fat_filament_pool  = [];
end

if (~isfield(output_feature, 'density_filament'))
    output_feature.density_filament  = [];
end

if (~isfield(output_feature, 'scrabled_density_filament'))
    output_feature.scrabled_density_filament  = [];
end


if (~isfield(output_feature, 'mean_st_per_fat_filament_pool'))
    output_feature.mean_st_per_fat_filament_pool  = [];
end

if (~isfield(output_feature, 'mean_st_per_filament_pool'))
    output_feature.mean_st_per_filament_pool  = [];
end

if (~isfield(output_feature, 'mean_intensity_per_fat_filament_pool'))
    output_feature.mean_intensity_per_fat_filament_pool  = [];
end

if (~isfield(output_feature, 'mean_intensity_per_filament_pool'))
    output_feature.mean_intensity_per_filament_pool  = [];
end


if (~isfield(output_feature, 'filament_mean_curvature'))
    output_feature.filament_mean_curvature  = [];
end

if (~isfield(output_feature, 'profileAllCell'))
    output_feature.profileAllCell  = [];
end

if (~isfield(output_feature, 'Centripetal_pixel'))
    
    output_feature.Centripetal_pixel  = [];
end



if (~isfield(output_feature, 'Cell_Mask'))
    output_feature.Cell_Mask  = [];
end

output_feature_reorganized = output_feature;

% arrange all the other feature to enable same format pooling

output_feature_reorganized.straightness_per_filament_pool = output_feature.straightness_per_filament_pool';
output_feature_reorganized.length_per_filament_pool = output_feature.length_per_filament_pool';
output_feature_reorganized.pixel_number_per_filament_pool = output_feature.pixel_number_per_filament_pool';

output_feature_reorganized.filament_mean_curvature = output_feature.filament_mean_curvature';

output_feature_reorganized.intensity_per_filament_pool = output_feature.intensity_per_filament_pool';
output_feature_reorganized.mean_intensity_per_filament_pool = output_feature.mean_intensity_per_filament_pool';
output_feature_reorganized.intensity_per_fat_filament_pool = output_feature.intensity_per_fat_filament_pool';
output_feature_reorganized.mean_intensity_per_fat_filament_pool = output_feature.mean_intensity_per_fat_filament_pool';

output_feature_reorganized.scale_per_filament_pool = output_feature.scale_per_filament_pool';

output_feature_reorganized.st_per_filament_pool = output_feature.st_per_filament_pool';
output_feature_reorganized.mean_st_per_filament_pool = output_feature.mean_st_per_filament_pool';
output_feature_reorganized.st_per_fat_filament_pool = output_feature.st_per_fat_filament_pool';
output_feature_reorganized.mean_st_per_fat_filament_pool = output_feature.mean_st_per_fat_filament_pool';

output_feature_reorganized.Centripetal_fila = output_feature.Centripetal_fila';

if(~isempty(output_feature.Cell_Mask) && ~isempty(output_feature.density_filament))
    output_feature_reorganized.density_filament = output_feature.density_filament(output_feature.Cell_Mask>0);
end

if(~isempty(output_feature.Cell_Mask) && ~isempty(output_feature.scrabled_density_filament))
    output_feature_reorganized.scrabled_density_filament = output_feature.scrabled_density_filament(output_feature.Cell_Mask>0);
end


