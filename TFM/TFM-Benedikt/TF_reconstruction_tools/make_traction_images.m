function make_traction_images(images_dir,TFM_results,TFM_settings, colorscale)
        %this program saves the traction data as fig and as tif. The
        %colorscale argument is optional and if given, it should be in in
        %the format: [lower limit, upper limit] e.g. [0,1000].
        
        result_dir = 'TFM_result_images';
        
        dir_struct = vertcat(dir(fullfile(images_dir,'*.tif*')));
        [sorted_dir_names,si] = sortrows({dir_struct.name}');
        stacklen = size(sorted_dir_names,1)
        if size(TFM_results,2) ~= stacklen
            disp('The number of images in the given directory is not equal to the number of traction datasets. Aborting!');
            return;
        end
        mkdir(result_dir);
        h1 = figure;
        h2 = figure;
        
        pack;
        for t = 1:stacklen
            bild = imread(fullfile(images_dir,sorted_dir_names{t}));
                       
            figure(h1); imagesc(bild); hold on; colormap gray, axis equal; 
            quiver(TFM_results(t).pos(:,1)+TFM_results(t).vec(:,1),TFM_results(t).pos(:,2)+TFM_results(t).vec(:,2),TFM_results(t).traction(:,1),TFM_results(t).traction(:,2),3,'r');
            title('Full stress data with locations corrected for substrate displacement');
            saveas(h1,fullfile(result_dir,['traction_vec_',num2str(t),'.tif']),'tiffn'); 
            saveas(h1,fullfile(result_dir,['traction_vec_',num2str(t),'.fig']),'fig');hold off;
            
            [grid_mat,tmat, i_max, j_max] = interp_vec2grid(TFM_results(t).pos+TFM_results(t).vec, TFM_results(t).traction, TFM_settings.meshsize); 
            tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;
            
            figure(h2); colormap default; surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm),view(0,90), shading interp, axis equal;
            set(gca, 'DataAspectRatio', [1,1,10],'YDir','reverse');
            if nargin >3
                set(gca,'CLim',colorscale);
            end
            colorbar;
            saveas(h2,fullfile(result_dir,['traction_mag_',num2str(t),'.tif']),'tiffn');
            saveas(h2,fullfile(result_dir,['traction_mag_',num2str(t),'.fig']),'fig');
        end
        close(h1);close(h2);
        pack;