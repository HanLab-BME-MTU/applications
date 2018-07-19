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

for iCompleteFrame = 1 : nCompleteFrame
    %     current_seg = current_seg_cell{1,iFrame};
    iFrame = iCompleteFrame+FirstFrame-1;
    
    smoothed_current_mask = smoothed_mask_cell{1,iCompleteFrame};
    
    C_xy = regionprops(smoothed_current_mask,'Centroid');
    
    center_x(iCompleteFrame) = C_xy.Centroid(1);
    center_y(iCompleteFrame) = C_xy.Centroid(2);
    
    current_VIF_image = MD.channels_(VIF_channel). loadImage(iFrame);
    
        RG_framem1 = zeros(size(smoothed_current_mask,1),size(smoothed_current_mask,2),3);
        RG_framem1(:,:,1) = smoothed_current_mask;
        
        cell_vif_pool = [cell_vif_pool current_VIF_image(smoothed_current_mask>0)];
        
        if iCompleteFrame>1
            RG_framem1(:,:,2) = smoothed_mask_cell{1,iCompleteFrame-1};
            previous_VIF_image = MD.channels_(VIF_channel). loadImage(iFrame-1);
            
            red_new = RG_framem1(:,:,1)-RG_framem1(:,:,2);
            green_old = RG_framem1(:,:,2)-RG_framem1(:,:,1);
            
            % for protrustion
            try
                red_vif_t_pool = [red_vif_t_pool current_VIF_image(red_new>0)];
            end
            
            % for retraction
            try
                green_vif_tm1_pool=[green_vif_tm1_pool previous_VIF_image(green_old>0)];
            end
            
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
    
    h5 = figure(5);    
    
    subplot(121); imagesc(RG_framem1);axis image; axis off;
    title('Difference','FontSize',15);
    axis([min_x max_x min_y max_y]);
    
    subplot(122); imagesc((region_branch_label_RGB));
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
    
    % find the vif intensity information
    for iL = 1 : trackedBranches        
        vif_pixel_values = current_VIF_image(find(labelMask==iL));
        vif_mean_matrix(iCompleteFrame,iL) = mean(vif_pixel_values);
    end
    
    title('Branches','FontSize',15);
    axis([min_x max_x min_y max_y]);
    
    saveas(h5,[outputPath,'\tracked_skel_region_',num2str(iFrame),'.tif']);
    
end
BA_output.branch_vif_mean_intensity  = nanmean(vif_mean_matrix');

BA_output.protrusion_vif_mean_intensity =  mean(red_vif_t_pool);

BA_output.retraction_vif_mean_intensity =  mean(green_vif_tm1_pool);

BA_output.whole_cell_vif_mean_intensity =  mean(cell_vif_pool);
