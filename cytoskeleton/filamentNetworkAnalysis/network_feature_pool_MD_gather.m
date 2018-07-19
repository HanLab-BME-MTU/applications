function network_feature_MD_thisCh_wholepool = network_feature_pool_MD_gather(network_feature, iChannel, feature_flag)
% function network_feature_MD_thisCh_wholepool = put_pool_together(network_feature_MD_thisCh_wholepool, network_feature,iChannel)
% not a function, is a script

network_feature_MD_thisCh_wholepool = [];

for iFrame = 1 : size(network_feature,2)
    
    if(isempty(network_feature_MD_thisCh_wholepool))
        network_feature_MD_thisCh_wholepool.straightness_per_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.length_per_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.pixel_number_per_filament_pool  = [];
        
        network_feature_MD_thisCh_wholepool.density_filament  = [];
        network_feature_MD_thisCh_wholepool.scrabled_density_filament   = [];
        
        network_feature_MD_thisCh_wholepool.orientation_pixel_pool_display  = [];
        network_feature_MD_thisCh_wholepool.orientation_pixel_pool_display_center  = [];
        
        network_feature_MD_thisCh_wholepool.intensity_per_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.mean_intensity_per_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.intensity_per_fat_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.mean_intensity_per_fat_filament_pool  = [];
        
        network_feature_MD_thisCh_wholepool.scale_per_filament_pool  = [];
        
        network_feature_MD_thisCh_wholepool.st_per_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.mean_st_per_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.st_per_fat_filament_pool  = [];
        network_feature_MD_thisCh_wholepool.mean_st_per_fat_filament_pool  = [];
        
        network_feature_MD_thisCh_wholepool.curvature_per_pixel_pool  = [];
        network_feature_MD_thisCh_wholepool.filament_mean_curvature  = [];
        
        network_feature_MD_thisCh_wholepool.nmssum_ratio_pool  = [];
        network_feature_MD_thisCh_wholepool.nmsmean_ratio_pool  = [];
        network_feature_MD_thisCh_wholepool.intsum_ratio_pool  = [];
        network_feature_MD_thisCh_wholepool.intmean_ratio_pool  = [];
        network_feature_MD_thisCh_wholepool.filasum_ratio_pool  = [];
        network_feature_MD_thisCh_wholepool.filamean_ratio_pool  = [];
        
        network_feature_MD_thisCh_wholepool.nmssum_ratio_pool_autoDpc  = [];
        network_feature_MD_thisCh_wholepool.nmsmean_ratio_pool_autoDpc  = [];
        network_feature_MD_thisCh_wholepool.intsum_ratio_pool_autoDpc  = [];
        network_feature_MD_thisCh_wholepool.intmean_ratio_pool_autoDpc  = [];
        network_feature_MD_thisCh_wholepool.filasum_ratio_pool_autoDpc  = [];
        network_feature_MD_thisCh_wholepool.filamean_ratio_pool_autoDpc  = [];
        
        network_feature_MD_thisCh_wholepool.Centripetal_fila  = [];
        network_feature_MD_thisCh_wholepool.Centripetal_pixel  = [];
        network_feature_MD_thisCh_wholepool.number_of_nucleus  = [];
        network_feature_MD_thisCh_wholepool.filament_density_mean = [];
        network_feature_MD_thisCh_wholepool.scrabled_density_filament_mean = [];
        
        network_feature_MD_thisCh_wholepool.ST_seg_ratio = [];
        network_feature_MD_thisCh_wholepool.GM_seg_ratio = [];
                       
    end
    
    if(~isempty(network_feature{iChannel, iFrame}))
        if(feature_flag(1)==1)
            network_feature_MD_thisCh_wholepool.straightness_per_filament_pool  = [ network_feature_MD_thisCh_wholepool.straightness_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.straightness_per_filament_pool(:)];
        end
        if(feature_flag(2)==1)
            network_feature_MD_thisCh_wholepool.length_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.length_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.length_per_filament_pool(:)];
        end
        if(feature_flag(3)==1)
            network_feature_MD_thisCh_wholepool.pixel_number_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.pixel_number_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.pixel_number_per_filament_pool(:)];
        end
        if(feature_flag(4)==1)
            
            network_feature_MD_thisCh_wholepool.density_filament    = [ network_feature_MD_thisCh_wholepool.density_filament;  ...
                network_feature{iChannel, iFrame}.density_filament(~isnan(network_feature{iChannel, iFrame}.density_filament))];
        end
        if(feature_flag(5)==1)
            network_feature_MD_thisCh_wholepool.scrabled_density_filament     = [ network_feature_MD_thisCh_wholepool.scrabled_density_filament;  ...
                network_feature{iChannel, iFrame}.scrabled_density_filament(~isnan(network_feature{iChannel, iFrame}.scrabled_density_filament))];
        end
        if(feature_flag(6)==1)
            network_feature_MD_thisCh_wholepool.orientation_pixel_pool_display    = [ network_feature_MD_thisCh_wholepool.orientation_pixel_pool_display;  ...
                network_feature{iChannel, iFrame}.orientation_pixel_pool_display(:)];
        end
        if(feature_flag(7)==1)
            network_feature_MD_thisCh_wholepool.orientation_pixel_pool_display_center    = [ network_feature_MD_thisCh_wholepool.orientation_pixel_pool_display_center;  ...
                network_feature{iChannel, iFrame}.orientation_pixel_pool_display_center(:)];
        end
        if(feature_flag(8)==1)
            network_feature_MD_thisCh_wholepool.intensity_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.intensity_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.intensity_per_filament_pool(:)];
        end
        if(feature_flag(9)==1)
            network_feature_MD_thisCh_wholepool.mean_intensity_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.mean_intensity_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.mean_intensity_per_filament_pool(:)];
        end
        if(feature_flag(10)==1)
            network_feature_MD_thisCh_wholepool.intensity_per_fat_filament_pool    = [ network_feature_MD_thisCh_wholepool.intensity_per_fat_filament_pool;  ...
                network_feature{iChannel, iFrame}.intensity_per_fat_filament_pool(:)];
        end
        if(feature_flag(11)==1)
            network_feature_MD_thisCh_wholepool.mean_intensity_per_fat_filament_pool    = [ network_feature_MD_thisCh_wholepool.mean_intensity_per_fat_filament_pool;  ...
                network_feature{iChannel, iFrame}.mean_intensity_per_fat_filament_pool(:)];
        end
        if(feature_flag(12)==1)
            network_feature_MD_thisCh_wholepool.scale_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.scale_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.scale_per_filament_pool(:)];
        end
        if(feature_flag(13)==1)
            network_feature_MD_thisCh_wholepool.st_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.st_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.st_per_filament_pool(:)];
        end
        if(feature_flag(14)==1)
            network_feature_MD_thisCh_wholepool.mean_st_per_filament_pool    = [ network_feature_MD_thisCh_wholepool.mean_st_per_filament_pool;  ...
                network_feature{iChannel, iFrame}.mean_st_per_filament_pool(:)];
        end
        if(feature_flag(15)==1)
            network_feature_MD_thisCh_wholepool.st_per_fat_filament_pool    = [ network_feature_MD_thisCh_wholepool.st_per_fat_filament_pool;  ...
                network_feature{iChannel, iFrame}.st_per_fat_filament_pool(:)];
        end
        if(feature_flag(16)==1)
            network_feature_MD_thisCh_wholepool.mean_st_per_fat_filament_pool    = [ network_feature_MD_thisCh_wholepool.mean_st_per_fat_filament_pool;  ...
                network_feature{iChannel, iFrame}.mean_st_per_fat_filament_pool(:)];
        end
        if(feature_flag(17)==1)
            network_feature_MD_thisCh_wholepool.curvature_per_pixel_pool    = [ network_feature_MD_thisCh_wholepool.curvature_per_pixel_pool;  ...
                network_feature{iChannel, iFrame}.curvature_per_pixel_pool(:)];
        end
        if(feature_flag(18)==1)
            network_feature_MD_thisCh_wholepool.filament_mean_curvature    = [ network_feature_MD_thisCh_wholepool.filament_mean_curvature;  ...
                network_feature{iChannel, iFrame}.filament_mean_curvature(:)];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.nmssum_ratio_pool    = [ network_feature_MD_thisCh_wholepool.nmssum_ratio_pool;  ...
            %         network_feature{iChannel, iFrame}.nmssum_ratio_pool];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.nmsmean_ratio_pool    = [ network_feature_MD_thisCh_wholepool.nmsmean_ratio_pool;  ...
            %         network_feature{iChannel, iFrame}.nmsmean_ratio_pool];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.intsum_ratio_pool    = [ network_feature_MD_thisCh_wholepool.intsum_ratio_pool;  ...
            %         network_feature{iChannel, iFrame}.intsum_ratio_pool];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.intmean_ratio_pool    = [ network_feature_MD_thisCh_wholepool.intmean_ratio_pool;  ...
            %         network_feature{iChannel, iFrame}.intmean_ratio_pool];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.filasum_ratio_pool    = [ network_feature_MD_thisCh_wholepool.filasum_ratio_pool;  ...
            %         network_feature{iChannel, iFrame}.filasum_ratio_pool];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.filamean_ratio_pool    = [ network_feature_MD_thisCh_wholepool.filamean_ratio_pool;  ...
            %         network_feature{iChannel, iFrame}.filamean_ratio_pool];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.nmssum_ratio_pool_autoDpc    = [ network_feature_MD_thisCh_wholepool.nmssum_ratio_pool_autoDpc;  ...
            %         network_feature{iChannel, iFrame}.nmssum_ratio_pool_autoDpc];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.nmsmean_ratio_pool_autoDpc    = [ network_feature_MD_thisCh_wholepool.straightness_per_filament_pool;  ...
            %         network_feature{iChannel, iFrame}.straightness_per_filament_pool];
            %  end
            % if(feature_flag(1)==1)
            %    network_feature_MD_thisCh_wholepool.intsum_ratio_pool_autoDpc    = [ network_feature_MD_thisCh_wholepool.intsum_ratio_pool_autoDpc;  ...
            %         network_feature{iChannel, iFrame}.intsum_ratio_pool_autoDpc];
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.intmean_ratio_pool_autoDpc    = [ network_feature_MD_thisCh_wholepool.intmean_ratio_pool_autoDpc;  ...
            %         network_feature{iChannel, iFrame}.intmean_ratio_pool_autoDpc];
            %
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.filasum_ratio_pool_autoDpc    = [ network_feature_MD_thisCh_wholepool.filasum_ratio_pool_autoDpc;  ...
            %         network_feature{iChannel, iFrame}.filasum_ratio_pool_autoDpc];
            %
            % end
            % if(feature_flag(1)==1)
            %     network_feature_MD_thisCh_wholepool.filamean_ratio_pool_autoDpc    = [ network_feature_MD_thisCh_wholepool.filamean_ratio_pool_autoDpc;  ...
            %         network_feature{iChannel, iFrame}.filamean_ratio_pool_autoDpc];
            % end
            if(feature_flag(21)==1)
                network_feature_MD_thisCh_wholepool.Centripetal_fila    = [ network_feature_MD_thisCh_wholepool.Centripetal_fila;  ...
                    network_feature{iChannel, iFrame}.Centripetal_fila(:)];
            end
            if(feature_flag(22)==1)
                network_feature_MD_thisCh_wholepool.Centripetal_pixel    = [ network_feature_MD_thisCh_wholepool.Centripetal_pixel;  ...
                    network_feature{iChannel, iFrame}.Centripetal_pixel(:)];
            end
                        
            
%             if(feature_flag(23)==1)
%                 network_feature_MD_thisCh_wholepool.Cell_Mask    =  [network_feature_MD_thisCh_wholepool.Cell_Mask;  ...
%                     network_feature{iChannel, iFrame}.Cell_Mask(:)];
%             end
            if(feature_flag(24)==1)
                network_feature_MD_thisCh_wholepool.number_of_nucleus    =  [network_feature_MD_thisCh_wholepool.number_of_nucleus;  ...
                    network_feature{iChannel, iFrame}.number_of_nucleus(:)];                
            end
            
            if(feature_flag(25)==1)
                network_feature_MD_thisCh_wholepool.filament_density_mean    =  [network_feature_MD_thisCh_wholepool.filament_density_mean;  ...
                    network_feature{iChannel, iFrame}.filament_density_mean(:)];
            end
            if(feature_flag(26)==1)
                network_feature_MD_thisCh_wholepool.scrabled_density_filament_mean    =  [network_feature_MD_thisCh_wholepool.scrabled_density_filament_mean;  ...
                    network_feature{iChannel, iFrame}.scrabled_density_filament_mean(:)];
            end
            
            if(feature_flag(27)==1)
                network_feature_MD_thisCh_wholepool.ST_seg_ratio    =  [network_feature_MD_thisCh_wholepool.ST_seg_ratio;  ...
                    network_feature{iChannel, iFrame}.ST_seg_ratio(:)];
            end
            if(feature_flag(28)==1)
                network_feature_MD_thisCh_wholepool.GM_seg_ratio    =  [network_feature_MD_thisCh_wholepool.GM_seg_ratio;  ...
                    network_feature{iChannel, iFrame}.GM_seg_ratio(:)];
            end
            
        end
    end
end
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
