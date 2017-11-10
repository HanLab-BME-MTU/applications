function ChMP_feature = extract_mean_percentiles_perwell_cell(NA_feature_thisMD, iChannel,iW, feature_index)

ChMP_feature = nan(18,8);

if(~isempty(NA_feature_thisMD{1}))
    for iF = 1:18
        this_feature=[];
        if(feature_index(iF)>0)
            for iFrame = (iW-1)*4+1 : (iW)*4                
                switch iF
                    case 1
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.straightness_per_filament_pool];
                    case 2
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.length_per_filament_pool];
                    case 3
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.pixel_number_per_filament_pool];
                    case 4
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.density_filament(NA_feature_thisMD{iChannel,iFrame}.Cell_Mask>0)'];
                    case 5
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.scrabled_density_filament(NA_feature_thisMD{iChannel,iFrame}.Cell_Mask>0)'];
                    case 6
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.orientation_pixel_pool_display'];
                    case 7
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.orientation_pixel_pool_display_center'];
                    case 8
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.intensity_per_filament_pool];
                    case 9
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.mean_intensity_per_filament_pool];
                    case 10
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.intensity_per_fat_filament_pool];
                    case 11
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.mean_intensity_per_fat_filament_pool];
                    case 12
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.scale_per_filament_pool];
                    case 13
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.st_per_filament_pool];
                    case 14
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.mean_st_per_filament_pool];
                    case 15
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.st_per_fat_filament_pool];
                    case 16
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.mean_st_per_fat_filament_pool];
                    case 17
                        try
                            this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.filament_mean_curvature];
                        end
                    case 18
                        try
                            this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.curvature_per_pixel_pool'];
                        end
                    otherwise
                        this_feature = [ this_feature NA_feature_thisMD{iChannel,iFrame}.straightness_per_filament_pool];
                end
                
                MP = mean_percentiles(this_feature);
                ChMP_feature(iF,1:8)=  MP(:);                
            end
        end
    end
    
end
function MP = mean_percentiles(this_feature)
this_feature = this_feature(:);
MP = nan(1,8);
% if there is nothing to get percentile, return nan
if( isempty(this_feature) || sum(isnan(this_feature))== numel(this_feature))
    return;
end

MP(1) = nanmean(this_feature(:));
MP(2) = prctile(this_feature(:),0);
MP(3) = prctile(this_feature(:),2);
MP(4) = prctile(this_feature(:),25);
MP(5) = prctile(this_feature(:),50);
MP(6) = prctile(this_feature(:),75);
MP(7) = prctile(this_feature(:),98);
MP(8) = prctile(this_feature(:),100);

