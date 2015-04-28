function output_feature_reorganized = arrange_network_features(output_feature_reorganized)

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


% arrange all the other feature to enable same format pooling

output_feature_reorganized.straightness_per_filament_pool = output_feature_reorganized.straightness_per_filament_pool';
output_feature_reorganized.length_per_filament_pool = output_feature_reorganized.length_per_filament_pool';
output_feature_reorganized.pixel_number_per_filament_pool = output_feature_reorganized.pixel_number_per_filament_pool';

output_feature_reorganized.filament_mean_curvature = output_feature_reorganized.filament_mean_curvature';

output_feature_reorganized.intensity_per_filament_pool = output_feature_reorganized.intensity_per_filament_pool';
output_feature_reorganized.mean_intensity_per_filament_pool = output_feature_reorganized.mean_intensity_per_filament_pool';
output_feature_reorganized.intensity_per_fat_filament_pool = output_feature_reorganized.intensity_per_fat_filament_pool';
output_feature_reorganized.mean_intensity_per_fat_filament_pool = output_feature_reorganized.mean_intensity_per_fat_filament_pool';

output_feature_reorganized.scale_per_filament_pool = output_feature_reorganized.scale_per_filament_pool';

output_feature_reorganized.st_per_filament_pool = output_feature_reorganized.st_per_filament_pool';
output_feature_reorganized.mean_st_per_filament_pool = output_feature_reorganized.mean_st_per_filament_pool';
output_feature_reorganized.st_per_fat_filament_pool = output_feature_reorganized.st_per_fat_filament_pool';
output_feature_reorganized.mean_st_per_fat_filament_pool = output_feature_reorganized.mean_st_per_fat_filament_pool';

output_feature_reorganized.Centripetal_fila = output_feature_reorganized.Centripetal_fila';

output_feature_reorganized.density_filament = output_feature_reorganized.density_filament(output_feature_reorganized.Cell_Mask>0);
output_feature_reorganized.scrabled_density_filament = output_feature_reorganized.scrabled_density_filament(output_feature_reorganized.Cell_Mask>0);

