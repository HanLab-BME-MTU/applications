% script for taking the 25%, 50% and 75% for each feature distriution
if(feature_flag(1)>0)
    quarters_straightness_per_filament_pool_all = [quarters_straightness_per_filament_pool_all;
        prctile(output_feature.straightness_per_filament_pool,[25 50 75])];
end

if(feature_flag(2)>0)
    quarters_length_per_filament_pool_all = [quarters_length_per_filament_pool_all;
        prctile(output_feature.length_per_filament_pool,[25 50 75])];
end

if(feature_flag(3)>0)
    quarters_pixel_number_per_filament_pool_all = [ quarters_pixel_number_per_filament_pool_all;
        prctile(output_feature.pixel_number_per_filament_pool,[25 50 75])];
end


if(feature_flag(4)>0 ||  feature_flag(5)>0 )
    % get masked density information
    density_masked = output_feature.density_filament(output_feature.Cell_Mask>0);
    quarters_density_filament_all = [ quarters_density_filament_all;
        prctile(density_masked,[25 50 75])];
    
    try
        scrabled_density_masked = output_feature.scrabled_density_filament(output_feature.Cell_Mask>0);
        
        quarters_scrabled_density_filament_all = [ quarters_scrabled_density_filament_all;
            prctile(scrabled_density_masked,[25 50 75])];
    end
end

if(feature_flag(6)>0 || feature_flag(7)>0)
    quarters_orientation_pixel_pool_display_all = [ quarters_orientation_pixel_pool_display_all;
        prctile(output_feature.orientation_pixel_pool_display,[25 50 75])];
    
    try
        quarters_orientation_pixel_pool_display_center_all = [ quarters_orientation_pixel_pool_display_center_all;
            prctile(output_feature.orientation_pixel_pool_display_center,[25 50 75])];
    end
end

if(feature_flag(8)>0)
    quarters_intensity_per_filament_pool_all = [ quarters_intensity_per_filament_pool_all;
        prctile(output_feature.intensity_per_filament_pool,[25 50 75])];    
end

if(feature_flag(9)>0)
    quarters_mean_intensity_per_filament_pool_all = [ quarters_mean_intensity_per_filament_pool_all;
        prctile(output_feature.mean_intensity_per_filament_pool,[25 50 75])];
end

if(feature_flag(10)>0)
    quarters_intensity_per_fat_filament_pool_all = [ quarters_intensity_per_fat_filament_pool_all;
        prctile(output_feature.intensity_per_fat_filament_pool,[25 50 75])];    
end

if(feature_flag(11)>0)
    quarters_mean_intensity_per_fat_filament_pool_all = [ quarters_mean_intensity_per_fat_filament_pool_all;
        prctile(output_feature.mean_intensity_per_fat_filament_pool,[25 50 75])];
end

if(feature_flag(12)>0)
    quarters_scale_per_filament_pool_all = [ quarters_scale_per_filament_pool_all;
        prctile(output_feature.scale_per_filament_pool,[25 50 75])];    
end

if(feature_flag(13)>0)
    quarters_st_per_filament_pool_all = [ quarters_st_per_filament_pool_all;
        prctile(output_feature.st_per_filament_pool,[25 50 75])];    
end

if(feature_flag(14)>0)
    quarters_mean_st_per_filament_pool_all = [ quarters_mean_st_per_filament_pool_all;
        prctile(output_feature.mean_st_per_filament_pool,[25 50 75])];    
end

if(feature_flag(15)>0)
    quarters_st_per_fat_filament_pool_all = [ quarters_st_per_fat_filament_pool_all;
        prctile(output_feature.st_per_fat_filament_pool,[25 50 75])];    
end

if(feature_flag(16)>0)
    quarters_mean_st_per_fat_filament_pooll_all = [ quarters_mean_st_per_fat_filament_pooll_all;
        prctile(output_feature.mean_st_per_fat_filament_pool,[25 50 75])];
end
