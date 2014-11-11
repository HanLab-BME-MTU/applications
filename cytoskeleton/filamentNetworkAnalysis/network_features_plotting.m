function network_features_plotting(output_feature, figure_flag, save_everything_flag,feature_flag,...
    im_name, outdir,iChannel,iFrame)
% function to plot network features
% Input: output_feature:   network feautre
%        figure_flag:      to plot to not
%        save_everything_flag: to save or not


% Plot only if the user requested

if(figure_flag>0)
    
    if(iscell(im_name))
        im_name = im_name{1};
    end
    
    % current_img = imread(im_name);
    
    im_name(im_name=='_')='-';
    im_name_title{1} = im_name;

    if(feature_flag(3)>0)
        % length distribution, by how many pixels
        h3 =  figure(3);
        
        [h,bin] = hist(output_feature.pixel_number_per_filament_pool,0:20:1000);
        h = h/length(output_feature.pixel_number_per_filament_pool);
        bar(bin,h);
        axis([-10 310 0 0.3]);
        
        title([im_name_title,' Pixels Number Distribution']);
        
        if(save_everything_flag>0)
            saveas(h3, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h3, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
    if(feature_flag(2)>0)
        
        h2 =  figure(2);
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
        h1 =  figure(1);
        [h,bin] = hist(output_feature.straightness_per_filament_pool,0:0.02:1);
        h = h/length(output_feature.straightness_per_filament_pool);
        bar(bin,h);
        axis([0.69 0.97 0 0.2]);
        
        title([im_name_title,' Straightness']);
        
        if(save_everything_flag>0)
            saveas(h1, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h1, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
        
        h5 =  figure(5);
        h = boxplot(output_feature.straightness_per_filament_pool);
        set(h(7,:),'Visible','off');
        axis([0 2 0.6 1.01]);
        
        title([im_name_title,' Straightness']);
        if(save_everything_flag>0)
            
            saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
       
    
    
    if(feature_flag(6)>0)
        
        
        h6 =  figure(6);
        rose(output_feature.orientation_pixel_pool_display);
        title([im_name_title,' Orientation of Filaments']);
        if(save_everything_flag>0)
            saveas(h6, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h6, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
        
        % the orientation in histogram
        h16 =  figure(16);
        
        [h,b] = hist(output_feature.orientation_pixel_pool_display,-pi/2:pi/18:pi/2);
        h=h(1:end);
        h = h/length(output_feature.orientation_pixel_pool_display);
        bin = (-pi/2:pi/18:pi/2) +pi/36;
        bar(bin,h);
        axis([-pi/2-pi/36 pi/2+pi/36 0 0.3]);
        
        title([im_name_title,' Orientation of Filaments']);
        if(save_everything_flag>0)
            saveas(h16, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h16, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
        end
    end
    
    if(feature_flag(7)>0)
        
        
        h7 =  figure(7);
        rose(output_feature.orientation_pixel_pool_display_center);
        title([im_name_title,' Orientation of Filaments']);
        
        if(save_everything_flag>0)
            saveas(h7, [outdir,filesep,'network_orientationrose_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h7, [outdir,filesep,'network_orientationrose_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
        
        % the orientation in histogram
        h17 =  figure(17);
        
        
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
            saveas(h17, [outdir,filesep,'network_orientationhist_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
            saveas(h17, [outdir,filesep,'network_orientationhist_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
        end
               
    end
end