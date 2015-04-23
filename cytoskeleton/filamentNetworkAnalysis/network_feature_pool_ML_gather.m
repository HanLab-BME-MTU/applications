function ML_pool_this_feature = network_feature_pool_ML_gather(NA_feature_thisMD, iChannel, feature_index)

if(~isempty(NA_feature_thisMD{iChannel}))
    for iF = 1:33
        this_feature=[];
        if(feature_index(iF)>0)
            switch iF
                case 1
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.straightness_per_filament_pool];
                case 2
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.length_per_filament_pool];
                case 3
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.pixel_number_per_filament_pool];
                case 4
                    %                        this_feature = [ this_feature NA_feature_thisMD{iChannel}.density_filament(NA_feature_thisMD{iChannel}.Cell_Mask>0)'];
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.density_filament];
                case 5
                    %                       this_feature = [ this_feature NA_feature_thisMD{iChannel}.scrabled_density_filament(NA_feature_thisMD{iChannel}.Cell_Mask>0)'];
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.scrabled_density_filament];
                case 6
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.orientation_pixel_pool_display];
                case 7
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.orientation_pixel_pool_display_center];
                case 8
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.intensity_per_filament_pool];
                case 9
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.mean_intensity_per_filament_pool];
                case 10
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.intensity_per_fat_filament_pool];
                case 11
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.mean_intensity_per_fat_filament_pool];
                case 12
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.scale_per_filament_pool];
                case 13
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.st_per_filament_pool];
                case 14
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.mean_st_per_filament_pool];
                case 15
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.st_per_fat_filament_pool];
                case 16
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.mean_st_per_fat_filament_pool];
                case 17
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.filament_mean_curvature];
                    end
                case 18
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.curvature_per_pixel_pool];
                    end
                case 19
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.nmssum_ratio_pool];
                    end
                case 20
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.nmsmean_ratio_pool];
                    end
                case 21
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.intsum_ratio_pool];
                    end
                case 22
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.intmean_ratio_pool];
                    end
                case 23
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.filasum_ratio_pool];
                    end
                case 24
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.filamean_ratio_pool];
                    end
                case 25
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.nmssum_ratio_pool_autoDpc];
                    end
                case 26
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.nmsmean_ratio_pool_autoDpc];
                    end
                case 27
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.intsum_ratio_pool_autoDpc];
                    end
                case 28
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.intmean_ratio_pool_autoDpc];
                    end
                case 29
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.filasum_ratio_pool_autoDpc];
                    end
                case 30
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.filamean_ratio_pool_autoDpc];
                    end
                case 31
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.Centripetal_fila];
                    end
                case 32
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.Centripetal_pixel];
                    end
                case 33
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.number_of_nucleus];
                    end
                case 34
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.filament_density_mean];
                    end
                    
                case 35
                    try
                        this_feature = [ this_feature; NA_feature_thisMD{iChannel}.scrabled_density_filament_mean];
                    end
                    
                otherwise
                    this_feature = [ this_feature; NA_feature_thisMD{iChannel}.straightness_per_filament_pool];
            end
            
        end
    end
    ML_pool_this_feature = this_feature;
else
    ML_pool_this_feature =[];
end