% function NA_feature_whole_thisMD = put_pool_together(NA_feature_whole_thisMD, NA_feature_thisMD,iChannel)
% not a function, is a script

if(isempty(NA_feature_whole_thisMD))
    NA_feature_whole_thisMD{iChannel} = [];
end

if(isempty(NA_feature_whole_thisMD{iChannel}))
    NA_feature_whole_thisMD{iChannel}.straightness_per_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.length_per_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.pixel_number_per_filament_pool  = [];
    
    NA_feature_whole_thisMD{iChannel}.density_filament  = [];
    NA_feature_whole_thisMD{iChannel}.scrabled_density_filament   = [];

    NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display  = [];
    NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display_center  = [];

    NA_feature_whole_thisMD{iChannel}.intensity_per_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.mean_intensity_per_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.intensity_per_fat_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.mean_intensity_per_fat_filament_pool  = [];
    
    NA_feature_whole_thisMD{iChannel}.scale_per_filament_pool  = [];
    
    NA_feature_whole_thisMD{iChannel}.st_per_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.mean_st_per_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.st_per_fat_filament_pool  = [];
    NA_feature_whole_thisMD{iChannel}.mean_st_per_fat_filament_pool  = [];

    NA_feature_whole_thisMD{iChannel}.curvature_per_pixel_pool  = [];
    NA_feature_whole_thisMD{iChannel}.filament_mean_curvature  = [];
    
    NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool  = [];
    NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool  = [];
    NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool  = [];
    NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool  = [];
    NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool  = [];
    NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool  = [];
    
    NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool_autoDpc  = [];
    NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool_autoDpc  = [];
    NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool_autoDpc  = [];
    NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool_autoDpc  = [];
    NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool_autoDpc  = [];
    NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool_autoDpc  = [];
    
    NA_feature_whole_thisMD{iChannel}.Centripetal_fila  = [];
    NA_feature_whole_thisMD{iChannel}.Centripetal_pixel  = [];
    NA_feature_whole_thisMD{iChannel}.number_of_nucleus  = [];
end


    
    
    NA_feature_whole_thisMD{iChannel}.straightness_per_filament_pool  = [ NA_feature_whole_thisMD{iChannel}.straightness_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.straightness_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.length_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.length_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.length_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.pixel_number_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.pixel_number_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.pixel_number_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.density_filament    = [ NA_feature_whole_thisMD{iChannel}.density_filament;  ...
        NA_feature_thisMD{iChannel, iFrame}.density_filament];
    NA_feature_whole_thisMD{iChannel}.scrabled_density_filament     = [ NA_feature_whole_thisMD{iChannel}.scrabled_density_filament;  ...
        NA_feature_thisMD{iChannel, iFrame}.scrabled_density_filament];
    NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display    = [ NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display;  ...
        NA_feature_thisMD{iChannel, iFrame}.orientation_pixel_pool_display];
    NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display_center    = [ NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display_center;  ...
        NA_feature_thisMD{iChannel, iFrame}.orientation_pixel_pool_display_center];
    NA_feature_whole_thisMD{iChannel}.intensity_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.intensity_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.intensity_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.mean_intensity_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.mean_intensity_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.mean_intensity_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.intensity_per_fat_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.intensity_per_fat_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.intensity_per_fat_filament_pool];
    NA_feature_whole_thisMD{iChannel}.mean_intensity_per_fat_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.mean_intensity_per_fat_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.mean_intensity_per_fat_filament_pool];
    NA_feature_whole_thisMD{iChannel}.scale_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.scale_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.scale_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.st_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.st_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.st_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.mean_st_per_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.mean_st_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.mean_st_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.st_per_fat_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.st_per_fat_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.st_per_fat_filament_pool];
    NA_feature_whole_thisMD{iChannel}.mean_st_per_fat_filament_pool    = [ NA_feature_whole_thisMD{iChannel}.mean_st_per_fat_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.mean_st_per_fat_filament_pool];
    NA_feature_whole_thisMD{iChannel}.curvature_per_pixel_pool    = [ NA_feature_whole_thisMD{iChannel}.curvature_per_pixel_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.curvature_per_pixel_pool];
    NA_feature_whole_thisMD{iChannel}.filament_mean_curvature    = [ NA_feature_whole_thisMD{iChannel}.filament_mean_curvature;  ...
        NA_feature_thisMD{iChannel, iFrame}.filament_mean_curvature];
    NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool    = [ NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.nmssum_ratio_pool];
    NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool    = [ NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.nmsmean_ratio_pool];
    NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool    = [ NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.intsum_ratio_pool];
    NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool    = [ NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.intmean_ratio_pool];
    NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool    = [ NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.filasum_ratio_pool];
    NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool    = [ NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.filamean_ratio_pool];
    NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool_autoDpc    = [ NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool_autoDpc;  ...
        NA_feature_thisMD{iChannel, iFrame}.nmssum_ratio_pool_autoDpc];
    NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool_autoDpc    = [ NA_feature_whole_thisMD{iChannel}.straightness_per_filament_pool;  ...
        NA_feature_thisMD{iChannel, iFrame}.straightness_per_filament_pool];
    NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool_autoDpc    = [ NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool_autoDpc;  ...
        NA_feature_thisMD{iChannel, iFrame}.intsum_ratio_pool_autoDpc];
    NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool_autoDpc    = [ NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool_autoDpc;  ...
        NA_feature_thisMD{iChannel, iFrame}.intmean_ratio_pool_autoDpc];
    
    NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool_autoDpc    = [ NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool_autoDpc;  ...
        NA_feature_thisMD{iChannel, iFrame}.filasum_ratio_pool_autoDpc];
    
    NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool_autoDpc    = [ NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool_autoDpc;  ...
        NA_feature_thisMD{iChannel, iFrame}.filamean_ratio_pool_autoDpc];
    NA_feature_whole_thisMD{iChannel}.Centripetal_fila    = [ NA_feature_whole_thisMD{iChannel}.Centripetal_fila;  ...
        NA_feature_thisMD{iChannel, iFrame}.Centripetal_fila];
    NA_feature_whole_thisMD{iChannel}.Centripetal_pixel    = [ NA_feature_whole_thisMD{iChannel}.Centripetal_pixel;  ...
        NA_feature_thisMD{iChannel, iFrame}.Centripetal_pixel];
    NA_feature_whole_thisMD{iChannel}.number_of_nucleus    =  [NA_feature_whole_thisMD{iChannel}.number_of_nucleus;  ...
        NA_feature_thisMD{iChannel, iFrame}.number_of_nucleus];
    
                                                       




% 
% textcell = {'Straightness of filament (for each filament)',... % 1
%     'Length (for each filament)',...                           % 2
%     'Pixel number of segmentation (for each filament)',...     % 3
%     'Filament density (for each valid pixel)',...                   % 4
%     'Scrabled filament density (for each valid pixel)',...          % 5
%     'Filament orientation (for each pixel)',...                % 6
%     'Filament orientation (for each pixel,centered)',...       % 7
%     'Filamet intensity (integrated for each filament)',...           % 8
%     'Filamet intensity (average for each filament)',...              % 9
%     'Filamet intensity (integrated for each dilated filament)',...   % 10
%     'Filamet intensity (average for each dilated filament)',...      % 11
%     'Scale of detected filaments (for each pixel)',...                   % 12
%     'ST Responce (integrated for each filament)',...            % 13
%     'ST Responce (average for each filament)',...               % 14
%     'ST Responce (integrated for each dilated filament)',...    % 15
%     'ST Responce (average for each dilated filament)',...       % 16
%     'Filament curvature (average for each filament)',...     % 17
%     'Filament curvature (for each pixel on filaments)',...   % 18
%     'NMS Sum Periphery-Center Ratio (for each sections for cell (Dpc40p))',...           % 19
%     'NMS Ave Periphery-Center Ratio (for each sections for cell(Dpc40p))',...            % 20
%     'Intensity Sum Periphery-Center Ratio (for each sections for cell(Dpc40p))',...      % 21
%     'Intensity Ave Periphery-Center Ratio (for each sections for cell(Dpc40p))',...      % 22
%     'Filament Pixel Sum Periphery-Center Ratio (for each sections for cell(Dpc40p))',... % 23
%     'Filament Density Periphery-Center Ratio (for each sections for cell(Dpc40p))',...   % 24
%     'NMS Sum Periphery-Center Ratio (for each sections for cell (AutoDpc))',...           % 25
%     'NMS Ave Periphery-Center Ratio (for each sections for cell(AutoDpc))',...            % 26
%     'Intensity Sum Periphery-Center Ratio (for each sections for cell(AutoDpc))',...      % 27
%     'Intensity Ave Periphery-Center Ratio (for each sections for cell(AutoDpc))',...      % 28
%     'Filament Pixel Sum Periphery-Center Ratio (for each sections for cell(AutoDpc))',... % 29
%     'Filament Density Periphery-Center Ratio (for each sections for cell(AutoDpc))',...   % 30
%     'Filament Centripetal angle (for each filament)',...   % 31
%     'Filament Centripetal angle (for each pixel)',...      % 32
%     'Number of Nucleus(per well)'...                    %33
%     };
% 

            
% %1
% straightness_per_filament_pool
% %2
% length_per_filament_pool: [1x1611 double]
% %3
% pixel_number_per_filament_pool: [1x1611 double]
% 
% %4
% density_filament: [2160x2560 double]
% %5
% scrabled_density_filament: [2160x2560 double]
% 
% %6
% orientation_pixel_pool_display: [43174x1 double]
% %7
% orientation_pixel_pool_display_center: [43174x1 double]
% 
% %8
% intensity_per_filament_pool: [1x1611 double]
% %9
% mean_intensity_per_filament_pool: [1x1611 double]
% %10
% intensity_per_fat_filament_pool: [1x1611 double]
% %11
% mean_intensity_per_fat_filament_pool: [1x1611 double]
% 
% %12
% scale_per_filament_pool: [1x1611 double]
% 
% %13
% st_per_filament_pool: [1x1611 double]
% %14
% mean_st_per_filament_pool: [1x1611 double]
% %15
% st_per_fat_filament_pool: [1x1611 double]
% %16
% mean_st_per_fat_filament_pool: [1x1611 double]
% 
% %17
% curvature_per_pixel_pool: [43174x1 double]
% %18
% filament_mean_curvature: [1x1611 double]
% 
% 
% 
% %19
% output_feature_reorganized.nmssum_ratio_pool = output_feature.profileAllCell.nmssum_pool;
% %20
% output_feature_reorganized.nmsmean_ratio_pool = output_feature.profileAllCell.nmsmean_pool;
% %21
% output_feature_reorganized.intsum_ratio_pool = output_feature.profileAllCell.intmean_pool;
% %22
% output_feature_reorganized.intmean_ratio_pool = output_feature.profileAllCell.intmean_pool;
% %23
% output_feature_reorganized.filasum_ratio_pool = output_feature.profileAllCell.filasum_pool;
% %24
% output_feature_reorganized.filamean_ratio_pool = output_feature.profileAllCell.filamean_pool;
% 
% 
% %25
% output_feature_reorganized.nmssum_ratio_pool_autoDpc   = profileStnmsSumPerpCenterRatio(:);
% %26
% output_feature_reorganized.nmsmean_ratio_pool_autoDpc  = profileStnmsMeanPerpCenterRatio(:);
% %27
% output_feature_reorganized.intsum_ratio_pool_autoDpc   = profileIntSumPerpCenterRatio(:);
% %28
% output_feature_reorganized.intmean_ratio_pool_autoDpc  = profileIntMeanPerpCenterRatio(:);
% %29
% output_feature_reorganized.filasum_ratio_pool_autoDpc  = profileFilaSumPerpCenterRatio(:);
% %30
% output_feature_reorganized.filamean_ratio_pool_autoDpc = profileFilaMeanPerpCenterRatio(:);
% 
% 
% %31
% Centripetal_fila: [1x1611 double]
% %32
% Centripetal_pixel: [41563x1 double]
% 
% %33
% number_of_nucleus: 30
% 
% 
% profileCell: {1x30 cell}
% profileAllCell: [1x1 struct]
% Cell_Mask: [2160x2560 double]
