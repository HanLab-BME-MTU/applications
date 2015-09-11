function network_features_plotting(output_feature, figure_flag, save_everything_flag,feature_flag,vimscreen_flag,...
    im_name, outdir,iChannel,iFrame,set_visible)
% function to plot network features
% Input: output_feature:   network feautre
%        figure_flag:      to plot to not
%        save_everything_flag: to save or not


% Plot only if the user requested
if isnumeric(set_visible)
    set_visible = 'on';
end
if(figure_flag>0)
    display(' --- Feature plotting');
    
    if(iscell(im_name))
        im_name = im_name{1};
    end
    
    % current_img = imread(im_name);
    
    im_name(im_name=='_')='-';
    im_name_title{1} = im_name;
    
    if(feature_flag(3)>0)
        % length distribution, by how many pixels
        
        h3 = figure(3); set(h3,'Visible',set_visible); hold off;
        
        [h,bin] = hist(output_feature.pixel_number_per_filament_pool,0:5:1000);
        h = h/length(output_feature.pixel_number_per_filament_pool);
        bar(bin,h);
        axis([-10 150 0 0.3]);
        
        title([im_name_title,' Pixels Number Distribution']);
        
        if(save_everything_flag>0)
            saveas(h3, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h3, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
    if(feature_flag(2)>0)
        
        
        h2 = figure(2); set(h2,'Visible',set_visible); hold off;
        
        [h,bin] = hist(output_feature.length_per_filament_pool,0:20:1000);
        h = h/length(output_feature.length_per_filament_pool);
        
        bar(bin,h);
        axis([-10 310 0 0.3]);
        
        title([im_name_title,' Length Distribution']);
        if(save_everything_flag>0)
            saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
    
    if(feature_flag(1)>0)
        
        % the straightness distribution
        
        h1 = figure(1); set(h1,'Visible',set_visible); hold off;
        
        [h,bin] = hist(output_feature.straightness_per_filament_pool,0:0.02:1);
        h = h/length(output_feature.straightness_per_filament_pool);
        bar(bin,h);
        axis([0.69 0.97 0 0.2]);
        
        title([im_name_title,' Straightness']);
        
        if(save_everything_flag>0)
            saveas(h1, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h1, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
        
        
        h5 = figure(5); set(h5,'Visible',set_visible); hold off;
        
        h = boxplot(output_feature.straightness_per_filament_pool);        
        try
            set(h(7,:),'Visible','off');
        end
        
        axis([0 2 0.6 1.01]);
        
        title([im_name_title,' Straightness']);
        if(save_everything_flag>0)
            
            saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
    if(feature_flag(6)>0)
        
        
        h6 = figure(6); set(h6,'Visible',set_visible); hold off;
        if(~isempty(output_feature.orientation_pixel_pool_display))
            rose(output_feature.orientation_pixel_pool_display);
        end
        
        title([im_name_title,' Orientation of Filaments']);
        if(save_everything_flag>0)
            saveas(h6, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h6, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
        
        % the orientation in histogram
        h106 = figure(106); set(h106,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.orientation_pixel_pool_display,-pi/2:pi/18:pi/2);
        h=h(1:end);
        h = h/length(output_feature.orientation_pixel_pool_display);
        bin = (-pi/2:pi/18:pi/2) +pi/36;
        bar(bin,h);
        axis([-pi/2-pi/36 pi/2+pi/36 0 0.3]);
        
        title([im_name_title,' Orientation of Filaments']);
        if(save_everything_flag>0)
            saveas(h106, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h106, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
    if(feature_flag(7)>0)
        
        
        h7 = figure(7); set(h7,'Visible',set_visible); hold off;
        if(~isempty(output_feature.orientation_pixel_pool_display_center))            
            rose(output_feature.orientation_pixel_pool_display_center);
        end
        
        title([im_name_title,' Orientation of Filaments']);
        
        if(save_everything_flag>0)
            saveas(h7, [outdir,filesep,'network_orientationrose_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h7, [outdir,filesep,'network_orientationrose_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
        
        % the orientation in histogram
        
        h107 = figure(107); set(h107,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.orientation_pixel_pool_display,-pi/2:pi/18:pi/2);
        h=h(1:end);
        h = h/length(output_feature.orientation_pixel_pool_display);
        
        ind_max_h = find(h==max(h));
        ind_max_h = ind_max_h(1);
        mode_bin = b(ind_max_h)-pi/36;
        mode_bin = mod(mode_bin,pi);
        
        [h,b] = hist(output_feature.orientation_pixel_pool_display_center,-pi/2+mode_bin:pi/18:pi/2+mode_bin);
        h=h(1:end);
        h = h/length(output_feature.orientation_pixel_pool_display_center);
        bin = (-pi/2:pi/18:pi/2)+mode_bin+pi/36;
        bar(bin,h);
        axis([-pi/2-pi/36+mode_bin pi/2+pi/36+mode_bin 0 0.3]);
        
        title([im_name_title,' Orientation of Filaments']);
        if(save_everything_flag>0)
            saveas(h107, [outdir,filesep,'network_orientationhist_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h107, [outdir,filesep,'network_orientationhist_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(8)>0)
        %intensity_per_filament_pool
        
        
        h8 = figure(8); set(h8,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.intensity_per_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.intensity_per_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Sum Intensity per fila']);
        if(save_everything_flag>0)
            saveas(h8, [outdir,filesep,'sum_intensity_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h8, [outdir,filesep,'sum_intensity_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    
    if(feature_flag(9)>0)
        %mean_intensity_per_filament_pool
        
        
        h9 = figure(9); set(h9,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.mean_intensity_per_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.mean_intensity_per_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Average Intensity per fila']);
        if(save_everything_flag>0)
            saveas(h9, [outdir,filesep,'mean_intensity_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h9, [outdir,filesep,'mean_intensity_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    
    if(feature_flag(10)>0)
        %intensity_per_fat_filament_pool
        
        
        h10 = figure(10); set(h10,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.intensity_per_fat_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.intensity_per_fat_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Sum Intensity per fat fila']);
        if(save_everything_flag>0)
            saveas(h10, [outdir,filesep,'sum_intensity_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h10, [outdir,filesep,'sum_intensity_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(11)>0)
        %mean_intensity_per_fat_filament_pool
        
        
        h11 = figure(11); set(h11,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.mean_intensity_per_fat_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.mean_intensity_per_fat_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Average Intensity per fat fila']);
        if(save_everything_flag>0)
            saveas(h11, [outdir,filesep,'mean_intensity_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h11, [outdir,filesep,'mean_intensity_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(12)>0)
        %scale_per_filament_pool
        
        
        h12 = figure(12); set(h12,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.scale_per_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.scale_per_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Average scale per fila']);
        if(save_everything_flag>0)
            saveas(h12, [outdir,filesep,'mean_scale_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h12, [outdir,filesep,'mean_scale_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    %%
    if(feature_flag(13)>0)
        %st_per_filament_pool
        
        
        h13 = figure(13); set(h13,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.st_per_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.st_per_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Sum ST per fila']);
        if(save_everything_flag>0)
            saveas(h13, [outdir,filesep,'sum_ST_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h13, [outdir,filesep,'sum_ST_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(14)>0)
        %st_per_filament_pool
        
        
        h14 = figure(14); set(h14,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.mean_st_per_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.mean_st_per_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Average ST per fila']);
        if(save_everything_flag>0)
            saveas(h14, [outdir,filesep,'mean_ST_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h14, [outdir,filesep,'mean_ST_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    if(feature_flag(15)>0)
        %st_per_fat_filament_pool
        
        
        h15 = figure(15); set(h15,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.st_per_fat_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.st_per_fat_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Sum ST per fat fila']);
        if(save_everything_flag>0)
            saveas(h15, [outdir,filesep,'sum_ST_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h15, [outdir,filesep,'sum_ST_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(16)>0)
        %mean_st_per_fat_filament_pool
        
        
        h16 = figure(16); set(h16,'Visible',set_visible); hold off;
        
        [h,b] = hist(output_feature.mean_st_per_fat_filament_pool);
        h=h(1:end);
        h = h/numel(output_feature.mean_st_per_fat_filament_pool);
        bar(b,h);
        
        title([im_name_title,' Average ST per fat fila']);
        if(save_everything_flag>0)
            saveas(h16, [outdir,filesep,'mean_ST_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h16, [outdir,filesep,'mean_ST_fatfila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    %%
    if(feature_flag(17)>0)
        %filament_mean_curvature
        h17 = figure(17); set(h17,'Visible',set_visible); hold off;
        
        [h,bin] = hist(output_feature.filament_mean_curvature);
        h = h/length(output_feature.filament_mean_curvature);
        bar(bin,h);
        
        title([im_name_title,' Ave Curvature per fila']);
        if(save_everything_flag>0)
            saveas(h17, [outdir,filesep,'ave_curvature_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h17, [outdir,filesep,'ave_curvature_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(18)>0)
        %curvature_per_pixel_pool
        h18 = figure(18); set(h18,'Visible',set_visible); hold off;
        
        [h,bin] = hist(abs(output_feature.curvature_per_pixel_pool),100);
        h = h/length(output_feature.curvature_per_pixel_pool);
        bar(bin,h);
        
        title([im_name_title,' Curvature per pixel']);
        if(save_everything_flag>0)
            saveas(h18, [outdir,filesep,'ave_curvature_pixel_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h18, [outdir,filesep,'ave_curvature_pixel_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    %%
    %     profileCell.nmsmean_pool
    
    if(feature_flag(20)>0 && vimscreen_flag>0)
        
        h20 = figure(20); set(h20,'Visible',set_visible); hold off;
        
        [h,bin] = hist((output_feature.profileAllCell.nmsmean_pool(:) ));
        h = h/length(output_feature.profileAllCell.nmsmean_pool(:) );
        bar(bin,h);
        title([im_name_title,' Average ST ratio(perp/center) ']);
        if(save_everything_flag>0)
            saveas(h20, [outdir,filesep,'mean_st_ratio_perp_center_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h20, [outdir,filesep,'mean_st_ratio_perp_center_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
        
        
        
        h120 = figure(120); set(h120,'Visible',set_visible); hold off;
        
        [h,bin] = hist((output_feature.profileAllCell.filamean_pool(:) ));
        h = h/length(output_feature.profileAllCell.filamean_pool(:) );
        bar(bin,h);
        title([im_name_title,' Average filament density ratio(perp/center) ']);
        if(save_everything_flag>0)
            saveas(h120, [outdir,filesep,'mean_st_ratio_perp_center_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h120, [outdir,filesep,'mean_st_ratio_perp_center_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    %%
    if(feature_flag(21)>0&& vimscreen_flag>0)
        %Centripetal_distribution
        
        h21 = figure(21); set(h21,'Visible',set_visible); hold off;
        
        [h,bin] = histc(output_feature.Centripetal_fila, 0:pi/24:pi/2);
        h(end-1) = h(end-1) +h(end);
        h = h(1:end-1);
        
        h = h/length(output_feature.Centripetal_fila);
        bar(pi/48:pi/24:pi/2-pi/48,h);
        
        title([im_name_title,' Centripetal per fila']);
        if(save_everything_flag>0)
            saveas(h21, [outdir,filesep,'centripetal_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h21, [outdir,filesep,'centripetal_fila_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(22)>0&&vimscreen_flag>0)
        %Centripetal_pixel_distribution
        
        h22 = figure(22); set(h22,'Visible',set_visible); hold off;
        
        [h,bin] = histc(output_feature.Centripetal_pixel, 0:pi/24:pi/2);
        h(end-1) = h(end-1) +h(end);
        h = h(1:end-1);
        h = h/length(output_feature.Centripetal_pixel);
        bar(pi/48:pi/24:pi/2-pi/48,h);
        
        title([im_name_title,' Centripetal per pixel']);
        if(save_everything_flag>0)
            saveas(h22, [outdir,filesep,'centripetal_pixel_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h22, [outdir,filesep,'centripetal_pixel_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(23)>0&&vimscreen_flag>0)
        %Centripetal_pixel_distribution
        
        h23 = figure(23); set(h23,'Visible',set_visible); hold off;
        
        hist(output_feature.number_of_nucleus);
        
        title([im_name_title,' Number of Cells']);
        if(save_everything_flag>0)
            saveas(h23, [outdir,filesep,'nucleus_number_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h23, [outdir,filesep,'nucleus_number_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(24)>0&&vimscreen_flag>0)
        %Centripetal_pixel_distribution
        
        h24 = figure(24); set(h24,'Visible',set_visible); hold off;
        
        hist(output_feature.filament_density_mean);
        
        title([im_name_title,' filament density mean']);
        
        if(save_everything_flag>0)
            saveas(h24, [outdir,filesep,'filament_density_mean_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h24, [outdir,filesep,'filament_density_mean_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
    if(feature_flag(25)>0&&vimscreen_flag>0)
        %Centripetal_pixel_distribution
        
        h25 = figure(25); set(h25,'Visible',set_visible); hold off;
        
        hist(output_feature.scrabled_density_filament_mean);
        
        title([im_name_title,'scrabled filament density mean']);
        
        if(save_everything_flag>0)
            saveas(h25, [outdir,filesep,'scrabled_filament_density_mean_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h25, [outdir,filesep,'scrabled_filament_density_mean_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
    end
    
end