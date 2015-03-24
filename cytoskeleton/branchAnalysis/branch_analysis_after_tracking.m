% for iCompleteFrame = 1 : nCompleteFrame
%     figure(6);imagesc(new_region_branch_label_cell{iCompleteFrame});
%     pause;
% end
%
% for iCompleteFrame = 1 : nCompleteFrame
%     imwrite(new_label_skel_cell{1,iCompleteFrame}, [outputPath,'\tracked_branch_',num2str(iCompleteFrame),'.tif']);
% end
trackedBranches=BA_output.branch_number_tracked;

red_vif_t_pool=[];
green_vif_tm1_pool=[];
cell_vif_pool = [];

cell_size_pool = [];
cell_vimtotal_pool = [];

cell_vif_seg_pool = [];
cell_vif_nms_pool = [];
  cell_vif_seg_total_pool=[];
  cell_vif_nms_total_pool=[];
  
% initialize empty pools
fila_branch_orientation_pool=[];
fila_trajectory_orientation_pool=[];
branch_trajectory_orientation_pool=[];


branch_filament_totallength_matrix=[];
branch_filament_meandensity_matrix=[];

branch_filament_totalnms_matrix=[];
branch_filament_meannms_matrix=[];

branch_orienation_perframe_allbranch = cell(1,1);

cell_vif_seg_total_array = [];
cell_vif_nms_total_array = [];

        

for iCompleteFrame = 1 :nCompleteFrame
%     current_seg = current_seg_cell{1,iCompleteFrame};
    iFrame = iCompleteFrame+FirstFrame-1;
    
    smoothed_current_mask = smoothed_mask_cell{1,iCompleteFrame};
    
    C_xy = regionprops(smoothed_current_mask,'Centroid');
    
    center_x(iCompleteFrame) = C_xy.Centroid(1);
    center_y(iCompleteFrame) = C_xy.Centroid(2);
    
    current_VIF_image = MD.channels_(VIF_channel). loadImage(iFrame);
    
    RG_framem1 = zeros(size(smoothed_current_mask,1),size(smoothed_current_mask,2),3);
    RG_framem1(:,:,1) = smoothed_current_mask;
    
    cell_vif_pool = [cell_vif_pool; current_VIF_image(smoothed_current_mask>0)];
    cell_size_pool = [cell_size_pool; sum(sum(smoothed_current_mask>0))];
    cell_vimtotal_pool = [cell_vimtotal_pool; sum(current_VIF_image(smoothed_current_mask>0))];
         
    
    if iCompleteFrame>1
        RG_framem1(:,:,2) = smoothed_mask_cell{1,iCompleteFrame-1};
        previous_VIF_image = MD.channels_(VIF_channel).loadImage(iFrame-1);
        
        red_new = RG_framem1(:,:,1)-RG_framem1(:,:,2);
        green_old = RG_framem1(:,:,2)-RG_framem1(:,:,1);
        
        % for protrustion
        %             try
        red_vif_t_pool = [red_vif_t_pool; current_VIF_image(red_new>0)];
        %             end
        
        % for retraction
        %             try
        green_vif_tm1_pool = [green_vif_tm1_pool; previous_VIF_image(green_old>0)];
        %             end
    else
        RG_framem1(:,:,2) = smoothed_current_mask;
    end
    
    blue_mask = current_mask*0;
    
    red_mask = RG_framem1(:,:,1);
    green_mask = RG_framem1(:,:,2);
    
    %%
    % Here is for color overlay
    %     red_mask(blue_mask>0)=0;
    %     green_mask(blue_mask>0)=0;
    %
    %     RG_framem1(:,:,1) = red_mask;
    %     RG_framem1(:,:,2) = green_mask;
    %
    %     RG_framem1(:,:,3) = blue_mask;
    
    %     subplot(223); imagesc(RG_framem1);axis image; axis off;
    %     title('Difference','FontSize',15);
    %     axis([min_x max_x min_y max_y]);
    
    
    labelMask = new_region_branch_label_cell{iCompleteFrame};
    
    region_branch_label = zeros(size(skel_no_branching));
    region_branch_label_R = region_branch_label;
    region_branch_label_G = region_branch_label;
    region_branch_label_B = region_branch_label;
    
    for iL = 1 : max(max(labelMask))
        %         region_branch_label==iL;
        region_branch_label_R(find(labelMask==iL))=(0.3+color_array(iL,1))/(1.3);
        region_branch_label_G(find(labelMask==iL))=(0.3+color_array(iL,2))/(1.3);
        region_branch_label_B(find(labelMask==iL))=(0.3+color_array(iL,3))/(1.3);
    end
    
    
    region_branch_label(smoothed_current_mask==0)=0;
    region_branch_label_R(smoothed_current_mask==0)=0;
    region_branch_label_G(smoothed_current_mask==0)=0;
    region_branch_label_B(smoothed_current_mask==0)=0;
    
    region_branch_label(smoothed_current_mask>0 & labelMask==0)=1000;
    region_branch_label_R(smoothed_current_mask>0 & labelMask==0)=0.5;
    region_branch_label_G(smoothed_current_mask>0 & labelMask==0)=0.5;
    region_branch_label_B(smoothed_current_mask>0 & labelMask==0)=0.5;
    
    region_branch_label_RGB = zeros( size(labelMask,1),size(labelMask,2),3);
    region_branch_label_RGB(:,:,1) = region_branch_label_R;
    region_branch_label_RGB(:,:,2) = region_branch_label_G;
    region_branch_label_RGB(:,:,3) = region_branch_label_B;
    
    region_orientation = region_orientation_cell{iCompleteFrame};
    skel_seg = (new_label_skel_cell{iCompleteFrame})>0;
    branch_only_orienation = region_orientation(skel_seg>0);
    
    if (iCompleteFrame<nCompleteFrame)
        trajectory_angle_this_frame = trajectory_angle(iCompleteFrame);
        branch_trajectory_orientation_pool = ...
            [branch_trajectory_orientation_pool; ...
            branch_only_orienation-trajectory_angle_this_frame];
    end
    
    % if there is filament segmentation
    if(~isempty(current_seg_cell{1,iCompleteFrame}) && filament_stat_flag>0)
        current_seg = current_seg_cell{1,iCompleteFrame};
        orienation_map_filtered = orienation_map_filtered_cell{1,iCompleteFrame};
        nms_map = nms_cell{1,iCompleteFrame};
        nms_map = nms_map.*current_seg;
        
        cell_vif_seg_total_pool = [cell_vif_seg_total_pool;...
            sum(sum(current_seg(smoothed_current_mask>0)))];
        cell_vif_nms_total_pool = [cell_vif_nms_total_pool; ...
            sum(sum(nms_map(smoothed_current_mask>0)))];
        cell_vif_seg_total_array(iCompleteFrame) = ...
            sum(sum(current_seg(smoothed_current_mask>0)));
        cell_vif_nms_total_array(iCompleteFrame) = ...
            sum(sum(nms_map(smoothed_current_mask>0)));
        
        AA = (pi/2-orienation_map_filtered.*current_seg);
        % wrap around in -pi/2 to pi/2
        AA(AA<-pi/2)=AA(AA<-pi/2)+pi;
        AA(AA<-pi/2)=AA(AA<-pi/2)+pi;
        AA(AA<-pi/2)=AA(AA<-pi/2)+pi;
        AA(AA>pi/2)=AA(AA>pi/2)-pi;
        AA(AA>pi/2)=AA(AA>pi/2)-pi;
        AA(AA>pi/2)=AA(AA>pi/2)-pi;
        
        filament_orientation = AA(current_seg>0);
        branch_orienation = region_orientation(current_seg>0);
        
        
        branch_orienation_perframe_allbranch{iCompleteFrame} = branch_orienation;
        
        fila_branch_orientation_pool = ...
            [fila_branch_orientation_pool; ...
            filament_orientation-branch_orienation];
        
        if (iCompleteFrame<nCompleteFrame)
            fila_trajectory_orientation_pool = ...
                [fila_trajectory_orientation_pool; ...
                filament_orientation-trajectory_angle_this_frame];
        end
    end
    
    if(figure_flag>0 && ~isempty(current_seg_cell{1,iCompleteFrame}) )
        
        [seg_ind_y,seg_ind_x] = find(current_seg>0 & labelMask > 0);
        
        h5 = figure(5); hold off;
        
        subplot(221); hold off; imagesc(RG_framem1);axis image; axis off;
        hold on;
        axis([min_x max_x min_y max_y]);
        
        
        hold on;
        plot(center_x(1:iCompleteFrame),center_y(1:iCompleteFrame),'r');
        hold on;
        plot(smooth_center_x(1:iCompleteFrame),smooth_center_y(1:iCompleteFrame));
        plot(smooth_center_x(1:iCompleteFrame),smooth_center_y(1:iCompleteFrame),'.');
        
        plot(center_x(1:iCompleteFrame),center_y(1:iCompleteFrame),'m.');
        
        plot(smooth_center_x(1),smooth_center_y(1),'*');
        plot(center_x(1),center_y(1),'m*');
        plot(smooth_center_x(iCompleteFrame),smooth_center_y(iCompleteFrame),'o');
        plot(center_x(iCompleteFrame),center_y(iCompleteFrame),'mo');
        
        try
        title({['Cell center speed: ', num2str(BA_output.speed_marked_frames(iCompleteFrame),'%.2f'), ' pixel/frame'], ...
            ['Smoothed speed: ', num2str(BA_output.smoothed_speed_marked_frames(iCompleteFrame),'%.2f'), ' pixel/frame'],...
           },'FontSize',15);
        end
        
         h105 = figure(105); hold off;
        
        hold off; imagesc(RG_framem1);axis image; axis off;
        hold on;
        axis([min_x max_x min_y max_y]);
                
        hold on;
        plot(center_x(1:iCompleteFrame),center_y(1:iCompleteFrame),'r');
        hold on;
        plot(smooth_center_x(1:iCompleteFrame),smooth_center_y(1:iCompleteFrame));
        plot(smooth_center_x(1:iCompleteFrame),smooth_center_y(1:iCompleteFrame),'.');
        
        plot(center_x(1:iCompleteFrame),center_y(1:iCompleteFrame),'m.');
        
        plot(smooth_center_x(1),smooth_center_y(1),'*');
        plot(center_x(1),center_y(1),'m*');
        plot(smooth_center_x(iCompleteFrame),smooth_center_y(iCompleteFrame),'o');
        plot(center_x(iCompleteFrame),center_y(iCompleteFrame),'mo');
        
        try
        title({['Cell center speed: ', num2str(BA_output.speed_marked_frames(iCompleteFrame),'%.2f'), ' pixel/frame'], ...
            ['Smoothed speed: ', num2str(BA_output.smoothed_speed_marked_frames(iCompleteFrame),'%.2f'), ' pixel/frame'],...
            ['Branch number: ', num2str(numel(unique(labelMask(:)))-1,'%.0f')]},'FontSize',15);
        end
        
        saveas(h105,[outputPath,filesep,'titled_tracked_region_',num2str(iFrame),'.tif']);
   
        h5=figure(5);
        
        subplot(222); hold off;
        imagesc((region_branch_label_RGB));
        axis image;axis off;
        hold on;
                
        for iL = 1 : trackedBranches
            
            [indy,indx]= find(new_label_skel_cell{iCompleteFrame}==iL);
            
            % find skel pixel in region map
            skel_region_values = labelMask(find(new_label_skel_cell{iCompleteFrame}==iL));
            if(mean(double(skel_region_values>0))<0.8)
                plot( indx,indy,'.','color',[ 0.2 0.2 0.2]);
            else
                plot( indx,indy,'.','color',color_array(iL,1:3)');
            end
        end
        hold on; plot(seg_ind_x,seg_ind_y,'b.','MarkerSize',1);    
        title({['Branch number: ', num2str(numel(unique(labelMask(:)))-1,'%.0f')]},'FontSize',15);       
        axis([min_x max_x min_y max_y]);
        
    end
    
%     % delete to release some memory
%     current_seg_cell{1,iCompleteFrame}=[];
%     orienation_map_filtered_cell{1,iCompleteFrame}=[];
    
    % find the vif intensity information
    for iL = 1 : trackedBranches
        vif_pixel_values = current_VIF_image(find(labelMask==iL));
        vif_mean_matrix(iCompleteFrame,iL) = mean(vif_pixel_values);
        
        if(sum(sum(current_seg))>0  && filament_stat_flag>0)
            vif_pixel_seg = current_seg(find(labelMask==iL));
            branch_filament_totallength_matrix(iCompleteFrame,iL) = sum(vif_pixel_seg);
            branch_filament_meandensity_matrix(iCompleteFrame,iL) = sum(vif_pixel_seg)./(numel(find(labelMask==iL)));
            
            vif_pixel_nms = nms_map(find(labelMask==iL));
            branch_filament_totalnms_matrix(iCompleteFrame,iL) = sum(vif_pixel_nms);
            branch_filament_meannms_matrix(iCompleteFrame,iL) = sum(vif_pixel_nms)./(numel(find(labelMask==iL)));
        end
        branch_size_matrix(iCompleteFrame,iL) = numel(find(labelMask==iL));
    end
    
    if(figure_flag>0 && ~isempty(current_seg_cell{1,iCompleteFrame}) && filament_stat_flag>0)    
        
        subplot(223); hold off;
        orient_display = branch_orienation;
               
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        
        rose(orient_display);
        title('Branch orientations','FontSize',15);
        
              
        h5 = figure(5); 
        subplot(224); hold off;
       
        orient_display = filament_orientation;
               
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        
        rose(orient_display);        
        title('Filament Orientations','FontSize',15);
        saveas(h5,[outputPath,filesep,'tracked_skel_region_',num2str(iFrame),'.tif']);
        
        
         h205 = figure(205); hold off;        
       
        orient_display = filament_orientation - branch_orienation;
               
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        orient_display(orient_display>pi/2)=...
            orient_display(orient_display>pi/2)-pi;
        orient_display(orient_display<-pi/2)=...
            orient_display(orient_display<-pi/2)+pi;
        
        rose(orient_display);        
        title(['Filament Orientations wrt Branch Orientation, std: ', num2str(std(orient_display),'%.2f')],'FontSize',15);
        saveas(h205,[outputPath,filesep,'fila_orient_wrt_branch_',num2str(iFrame),'.tif']);
  
    end
    
end

% if filament stat is requested and available
if(~isempty(fila_branch_orientation_pool) && filament_stat_flag>0)
    % wrap the angles in -pi/2 to pi/2
    fila_branch_orientation_pool(fila_branch_orientation_pool>pi/2)=...
        fila_branch_orientation_pool(fila_branch_orientation_pool>pi/2)-1*pi;
    fila_branch_orientation_pool(fila_branch_orientation_pool>pi/2)=...
        fila_branch_orientation_pool(fila_branch_orientation_pool>pi/2)-1*pi;
    fila_branch_orientation_pool(fila_branch_orientation_pool<-pi/2)=...
        fila_branch_orientation_pool(fila_branch_orientation_pool<-pi/2)+1*pi;
    fila_branch_orientation_pool(fila_branch_orientation_pool<-pi/2)=...
        fila_branch_orientation_pool(fila_branch_orientation_pool<-pi/2)+1*pi;
    
    
    if(figure_flag>0)
        % filament vs cellmovement
        [h,bin]= hist(fila_branch_orientation_pool,-pi/2+pi/36:pi/(18):pi/2-pi/36);
        h = h./(sum(h))*100;
        h12 = figure(12); hold off;
        bar(bin, h);
        axis([-pi/2 pi/2 0 max(h)+1]);
        set(gca, 'xtick', -pi/2:pi/4:pi/2);
        set(gca, 'xticklabel', {'-pi/2','-pi/4','0','pi/4','pi/2'});
        title('Orientation difference between filament and branch');
        xlabel('Orientation Difference (unit: rad)');
        ylabel('Percentage(%)');
        saveas(h12,[outputPath,filesep,'fila_vs_branch_stat.tif']);
        
       
    end
end

if(~isempty(fila_trajectory_orientation_pool) && filament_stat_flag>0)
    fila_trajectory_orientation_pool(fila_trajectory_orientation_pool>pi/2)=...
        fila_trajectory_orientation_pool(fila_trajectory_orientation_pool>pi/2)-1*pi;
    fila_trajectory_orientation_pool(fila_trajectory_orientation_pool>pi/2)=...
        fila_trajectory_orientation_pool(fila_trajectory_orientation_pool>pi/2)-1*pi;
    fila_trajectory_orientation_pool(fila_trajectory_orientation_pool<-pi/2)=...
        fila_trajectory_orientation_pool(fila_trajectory_orientation_pool<-pi/2)+1*pi;
    fila_trajectory_orientation_pool(fila_trajectory_orientation_pool<-pi/2)=...
        fila_trajectory_orientation_pool(fila_trajectory_orientation_pool<-pi/2)+1*pi;
    
    if(figure_flag>0)
        % filament vs cellmovement
        [h,bin]= hist(fila_trajectory_orientation_pool,-pi/2+pi/36:pi/(18):pi/2-pi/36);
        h = h./(sum(h))*100;
        h13 = figure(13); hold off;
        bar(bin, h);
        axis([-pi/2 pi/2 0 max(h)+1]);
        set(gca, 'xtick', -pi/2:pi/4:pi/2);
        set(gca, 'xticklabel', {'-pi/2','-pi/4','0','pi/4','pi/2'});
        title('Orientation difference between filament and cell movement');
        xlabel('Orientation Difference (unit: rad)');
        ylabel('Percentage(%)');
        saveas(h13,[outputPath,filesep,'fila_vs_cellmove_stat.tif']);
    end
end

if(~isempty(branch_trajectory_orientation_pool) && filament_stat_flag>0)
    % wrap the angles in -pi/2 to pi/2
    branch_trajectory_orientation_pool(branch_trajectory_orientation_pool>pi/2)=...
        branch_trajectory_orientation_pool(branch_trajectory_orientation_pool>pi/2)-1*pi;
    branch_trajectory_orientation_pool(branch_trajectory_orientation_pool>pi/2)=...
        branch_trajectory_orientation_pool(branch_trajectory_orientation_pool>pi/2)-1*pi;
    
    branch_trajectory_orientation_pool(branch_trajectory_orientation_pool<-pi/2)=...
        branch_trajectory_orientation_pool(branch_trajectory_orientation_pool<-pi/2)+1*pi;
    branch_trajectory_orientation_pool(branch_trajectory_orientation_pool<-pi/2)=...
        branch_trajectory_orientation_pool(branch_trajectory_orientation_pool<-pi/2)+1*pi;
    
    
    if(figure_flag>0)
        %plot branch vs cell movement distribution
        [h,bin]= hist(branch_trajectory_orientation_pool,-pi/2+pi/36:pi/(18):pi/2-pi/36);
        h = h./(sum(h))*100;
        h14 = figure(14); hold off;
        bar(bin, h);
        axis([-pi/2 pi/2 0 max(h)+1]);
        set(gca, 'xtick', -pi/2:pi/4:pi/2);
        set(gca, 'xticklabel', {'-pi/2','-pi/4','0','pi/4','pi/2'});
        title('Orientation difference between branch orientation and cell movement');
        xlabel('Orientation Difference (unit: rad)');
        ylabel('Percentage(%)');
        saveas(h14,[outputPath,filesep,'branch_vs_cellmove_stat.tif']);
    end
end

BA_output.branch_vif_mean_intensity  = nanmean(vif_mean_matrix);

%%

if(sum(sum(current_seg))>0  && filament_stat_flag>0)
    
    branch_filament_totallength_matrix(branch_filament_totallength_matrix==0)=nan;
    branch_filament_totalnms_matrix(branch_filament_totalnms_matrix==0)=nan;
    
    BA_output.branch_seg_total = nanmean(branch_filament_totallength_matrix);
    BA_output.branch_seg_mean = nanmean(branch_filament_meandensity_matrix);
    BA_output.branch_nms_total = nanmean(branch_filament_totalnms_matrix);
    BA_output.branch_nms_mean = nanmean(branch_filament_meannms_matrix);
else
    BA_output.branch_seg_total = nan;
    BA_output.branch_seg_mean = nan;
    BA_output.branch_nms_total = nan;
    BA_output.branch_nms_mean = nan;    
end


%%

BA_output.branch_filament_totallength_matrix = branch_filament_totallength_matrix;
BA_output.branch_filament_meandensity_matrix = branch_filament_meandensity_matrix;
BA_output.branch_filament_totalnms_matrix = branch_filament_totalnms_matrix;
BA_output.branch_filament_meannms_matrix = branch_filament_meannms_matrix;
BA_output.branch_vif_mean_matrix = vif_mean_matrix;
BA_output.branch_size_matrix = branch_size_matrix;

BA_output.cell_size_array = cell_size_pool;
BA_output.cell_vif_seg_total_array = cell_vif_seg_total_array;
BA_output.cell_vif_nms_total_array = cell_vif_nms_total_array;
BA_output.cell_vimtotal_pool = cell_vimtotal_pool;
BA_output.branch_orienation_perframe_allbranch = branch_orienation_perframe_allbranch;
%%


BA_output.branch_mean_size  = nanmean(branch_size_matrix);

BA_output.protrusion_vif_mean_intensity =  mean(red_vif_t_pool);

BA_output.retraction_vif_mean_intensity =  mean(green_vif_tm1_pool);

BA_output.whole_cell_vif_mean_intensity =  mean(cell_vif_pool);

%%
if(sum(sum(current_seg))>0  && filament_stat_flag>0)
    
    BA_output.whole_cell_vim_seg_total = nanmean(cell_vif_seg_total_pool);
    BA_output.whole_cell_vim_seg_mean = nanmean(cell_vif_seg_total_pool./cell_size_pool);
    
    BA_output.whole_cell_vim_nms_total = nanmean(cell_vif_nms_total_pool);
    BA_output.whole_cell_vim_nms_mean = nanmean(cell_vif_nms_total_pool./cell_size_pool);
    
else
    BA_output.whole_cell_vim_seg_total =nan;
    BA_output.whole_cell_vim_seg_mean = nan;
    
    BA_output.whole_cell_vim_nms_total = nan;
    BA_output.whole_cell_vim_nms_mean = nan;
    
end
%%

%get a random sampling from all the vif_pool, by a fixed number of frames
%to be considered
BA_output.pool_all_vif_intensity = ...
    cell_vif_pool(randsample(numel(cell_vif_pool),...
    round(numel(cell_vif_pool)/max(1,(BA_output.cell_marked_frame_number/10)))));

BA_output.whole_cell_vim_totalamount_mean  = mean(cell_vimtotal_pool);

BA_output.whole_cell_size_mean  = mean(cell_size_pool);

% if not calculated, these are empty
BA_output.fila_branch_orientation_pool = fila_branch_orientation_pool;

BA_output.fila_branch_orientation_pool_std = std(fila_branch_orientation_pool);


BA_output.fila_trajectory_orientation_pool = fila_trajectory_orientation_pool;

BA_output.fila_trajectory_orientation_pool_std = std(fila_trajectory_orientation_pool);


BA_output.branch_trajectory_orientation_pool = branch_trajectory_orientation_pool;

BA_output.branch_cellmovement_std = std(branch_trajectory_orientation_pool);
