% for iCompleteFrame = 1 : nCompleteFrame
%     figure(6);imagesc(new_region_branch_label_cell{iCompleteFrame});
%     pause;
% end
%
% for iCompleteFrame = 1 : nCompleteFrame
%     imwrite(new_label_skel_cell{1,iCompleteFrame}, [outputPath,'\tracked_branch_',num2str(iCompleteFrame),'.tif']);
% end
trackedBranches=BA_output.branch_number_tracked;

radius = 40;
min_change = 100;
peri_dist = 60;

BA_output.curvature_map_pool = [];
BA_output.filament_alignment_map_pool = [];
BA_output.filament_density_pool = [];
BA_output.prot_or_retr_fila_pool = [];
BA_output.curvature_map_full_pool = [];
BA_output.filament_alignment_map_full_pool = [];
BA_output.filament_density_full_pool = [];
BA_output.prot_or_retr_fila_full_pool = [];
BA_output.scale_map_pool = [];

BA_output.last_curvature_map_pool = [];
BA_output.last_filament_alignment_map_pool = [];
BA_output.last_filament_density_pool = [];
BA_output.last_prot_or_retr_fila_pool = [];
BA_output.last_curvature_map_full_pool = [];
BA_output.last_filament_alignment_map_full_pool = [];
BA_output.last_filament_density_full_pool = [];
BA_output.last_prot_or_retr_fila_full_pool = [];
BA_output.last_scale_map_pool = [];

BA_output.int_map_pool = [];
BA_output.st_map_pool = [];
BA_output.int_map_full_pool = [];
BA_output.st_map_full_pool = [];

BA_output.last_int_map_pool = [];
BA_output.last_st_map_pool = [];
BA_output.last_int_map_full_pool = [];
BA_output.last_st_map_full_pool = [];

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



% for iCompleteFrame = 1 :nCompleteFrame-1
% %     current_seg = current_seg_cell{1,iCompleteFrame};
%     iFrame = iCompleteFrame+FirstFrame-1;
%     smoothed_current_mask = smoothed_mask_cell{1,iCompleteFrame};
%     current_VIF_image = MD.channels_(VIF_channel). loadImage(iFrame);
%
%     next_smoothed_mask = smoothed_mask_cell{1,iCompleteFrame+1};
%     next_VIF_image = MD.channels_(VIF_channel). loadImage(iFrame+1);
%
%     optical_flow_for_vif;
%     saveas(h101,[outputPath,filesep,'vif_flow_frame_',num2str(iFrame),'.tif']);
%
% end


this_frame_prot_retract_map_based_on_branch_cell = cell(1,nCompleteFrame);
last_frame_prot_retract_map_based_on_branch_cell = cell(1,nCompleteFrame);
filament_density_cell =  cell(1,nCompleteFrame);
curvature_map_cell =  cell(1,nCompleteFrame);
filament_alignment_map_cell =  cell(1,nCompleteFrame);

for iCompleteFrame = 1 :nCompleteFrame
    %     current_seg = current_seg_cell{1,iCompleteFrame};
    iFrame = iCompleteFrame+FirstFrame-1;
    
    % this frame, post protusion or retraction
    this_frame_prot_retract_map_based_on_branch = zeros(size(smoothed_mask_cell{1,1}));
    
    % previous frame, pre- protusion or retraction
    last_frame_prot_retract_map_based_on_branch = zeros(size(smoothed_mask_cell{1,1}));
    
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
        
        red_new_positive  = red_new>0;
        green_old_positive  = green_old>0;
        
        this_cell_mask = squeeze(RG_framem1(:,:,1));
        last_cell_mask = squeeze(RG_framem1(:,:,2));
        
        this_cell_mask_center = (bwdist(1-this_cell_mask))>peri_dist;
        last_cell_mask_center = (bwdist(1-last_cell_mask))>peri_dist;
        
        
        
        for iT = 1 : trackedBranches
            this_branch_this_frame = ...
                (new_region_branch_label_cell{iCompleteFrame}==iT);
            
            if(iCompleteFrame>1)
                this_branch_last_frame = ...
                    (new_region_branch_label_cell{iCompleteFrame-1}==iT);
            else
                this_branch_last_frame=[];
            end
            
            this_branch_min_change = min(min_change, ...
                sum(sum(this_branch_last_frame))/3);
            last_branch_min_change = min(min_change, ...
                sum(sum(this_branch_last_frame))/3);
            
            
            %              if(iFrame<=numel(new_region_branch_label_cell)-1)
            %                  this_branch_next_frame = ...
            %                      (new_region_branch_label_cell{iFrame+1}==iT);
            %              else
            %                  this_branch_next_frame=[];
            %              end
            
            this_frame_prot_retract_map_based_on_branch(this_branch_this_frame>0)=0.5;
            last_frame_prot_retract_map_based_on_branch(this_branch_last_frame>0)=0.5;
            
            if(sum(sum(this_branch_this_frame))>0)
                
                if(sum(sum(this_branch_last_frame))>0)
                    
                    % prot i-1 to i frame counted in frame i
                    % compared with
                    % retract from i-1 to i frame counted in frame i-1
                    
                    % and also need to be more than a number of pixels, if
                    % really just a few pixels, don't consider
                    
                    if(sum(sum(red_new_positive(this_branch_this_frame>0))) > ...
                            max(this_branch_min_change,sum(sum(green_old_positive(this_branch_last_frame>0)))) )
                        this_frame_prot_retract_map_based_on_branch(this_branch_this_frame>0)=1;
                        last_frame_prot_retract_map_based_on_branch(this_branch_last_frame>0)=1;
                    else
                        if(sum(sum(green_old_positive(this_branch_last_frame>0)))>last_branch_min_change)
                            this_frame_prot_retract_map_based_on_branch(this_branch_this_frame>0)=-1;
                            last_frame_prot_retract_map_based_on_branch(this_branch_last_frame>0)=-1;
                        end
                    end
                    
                else
                    % if there is no previous frame for this branch
                    if(sum(sum(red_new_positive(this_branch_this_frame>0))) > this_branch_min_change)
                        this_frame_prot_retract_map_based_on_branch(this_branch_this_frame>0)=1;
                    end
                end
            else
                last_frame_prot_retract_map_based_on_branch(this_branch_last_frame>0)=-1;
                
            end
        end
        
        % correction on direct pixel comparison
        
        this_frame_prot_retract_map_based_on_branch(red_new_positive>0) = 1;
        %         last_frame_prot_retract_map_based_on_branch(green_old_positive>0) = -1;
        
        % put the center away
        this_frame_prot_retract_map_based_on_branch(this_cell_mask_center>0) = 0.5;
        last_frame_prot_retract_map_based_on_branch(last_cell_mask_center>0) = 0.5;
        
        
        % plot for debugging
        
        h17=figure(17);
        subplot(121);
        imagesc(last_frame_prot_retract_map_based_on_branch);
        axis image; axis off; colormap(gray);caxis([-1 1]);
        subplot(122);
        imagesc(this_frame_prot_retract_map_based_on_branch);
        axis image; axis off; colormap(gray);caxis([-1 1]);
        saveas(h17, [outputPath,filesep,'pre_prot_retract_this_regions_',num2str(iFrame),'.tif']);
        saveas(h17, [outputPath,filesep,'pre_prot_retract_this_regions_',num2str(iFrame),'.fig']);
        
        
        % keep down the protusion and retraction definition
        this_frame_prot_retract_map_based_on_branch_cell{iCompleteFrame} = this_frame_prot_retract_map_based_on_branch;
        last_frame_prot_retract_map_based_on_branch_cell{iCompleteFrame} = last_frame_prot_retract_map_based_on_branch;
        
       
        
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
        this_frame_prot_retract_map_based_on_branch_cell{iCompleteFrame} = this_frame_prot_retract_map_based_on_branch;
        last_frame_prot_retract_map_based_on_branch_cell{iCompleteFrame} = last_frame_prot_retract_map_based_on_branch;
        
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
        
        % network analysis for the following 3 properties
        % 1. Straghtness using the curvalture definition
        % 2. Filament orientation local alignment
        % 3. Filament Compression -- density
        % 4. Liya added: vim intensity(5), ST(6)
        
        % 1, Straghtness, has to be calculated from model
        
        curvature_map = nan(size(smoothed_mask_cell{1,1}));
        
        current_model = current_model_cell{iCompleteFrame};
        for iFm = 1 : length(current_model)
            %             try
            x = (current_model{iFm}(:,1));
            y = (current_model{iFm}(:,2));
            
            line_smooth_H = fspecial('gaussian',11,1.5);
            
            line_i_x = (imfilter(x, line_smooth_H, 'replicate', 'same'));
            line_i_y = (imfilter(y, line_smooth_H, 'replicate', 'same'));
            
            Vertices = [line_i_x line_i_y];
            Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
            k=LineCurvature2D(Vertices,Lines);
            
            curvature_map(sub2ind(size(smoothed_mask_cell{1,1}),round(y),round(x))) = abs(k);
            %             end
        end
        
        curvature_map_cell{iCompleteFrame} = curvature_map;
        
        % 2. Orientation local alignment as circular std
        filament_orientation_map = AA;
        filament_orientation_map(current_seg==0)=nan;
        
        dilated_current_seg = current_seg;
        [ify_array, ifx_array] = find(current_seg>0);
        
        filament_alignment_map = nan(size(smoothed_mask_cell{1,1}));
        
        for if_ind = 1 : numel(ify_array)
            ifx = ifx_array(if_ind);
            ify = ify_array(if_ind);
            
            try
                filament_squre = filament_orientation_map( round(ify-radius/2): round(ify+radius/2), ...
                    round(ifx-radius/2): round(ifx+radius/2));
                filament_small_pool = filament_squre(isnan(filament_squre)==0);
                alignment_value = circ_std(filament_small_pool);
                filament_alignment_map(ify,ifx) = real(alignment_value);
            catch
                filament_alignment_map(ify,ifx) = nan;
            end
        end
        
        filament_alignment_map_cell{iCompleteFrame} = filament_alignment_map;
        
        % 3. Density
        current_seg_double = double(current_seg);
        
        filament_density = imfilter(current_seg_double, ones(radius, radius)/(radius*radius), 'same','replicate');
        
        filament_density_cell{iCompleteFrame} = filament_density;
        
        if(iCompleteFrame>1)
            this_nms = nms_cell{iCompleteFrame};
            this_scaleMap =  scaleMap_cell{1,iCompleteFrame};
            last_scaleMap =  scaleMap_cell{1,iCompleteFrame-1};
            
            last_VIF_image = MD.channels_(VIF_channel).loadImage(iFrame-1);
            last_nms = nms_cell{iCompleteFrame-1};
            last_MAX_st_res = MAX_st_res_cell{iCompleteFrame-1};
            this_MAX_st_res = MAX_st_res_cell{iCompleteFrame};
            
            curvature_map_pool = curvature_map(current_seg>0);
            filament_alignment_map_pool = filament_alignment_map(current_seg>0);
            filament_density_pool = filament_density(current_seg>0);
            scale_map_pool = this_scaleMap(current_seg>0);
            prot_or_retr_fila_pool = this_frame_prot_retract_map_based_on_branch(current_seg>0);
            
            
            curvature_map(isnan(curvature_map))=0;
            filament_alignment_map(isnan(filament_alignment_map))=0;
            
            
            curvature_map_full = imfilter(curvature_map,ones(radius, radius)/(radius*radius), 'same','replicate')./filament_density;
            filament_alignment_map_full = imfilter(filament_alignment_map,ones(radius, radius)/(radius*radius), 'same','replicate')./filament_density;
            
            curvature_map_full_pool = curvature_map_full(filament_density>0);
            filament_alignment_map_full_pool = filament_alignment_map_full(filament_density>0);
            filament_density_full_pool = filament_density(filament_density>0);
            prot_or_retr_fila_full_pool = this_frame_prot_retract_map_based_on_branch(filament_density>0);
            
            int_map_pool = current_VIF_image(current_seg>0);
            st_map_pool = this_nms(current_seg>0);
            
            int_map_full_pool = current_VIF_image(filament_density>0);
            st_map_full_pool = this_MAX_st_res(filament_density>0);
            
            BA_output.curvature_map_pool = [BA_output.curvature_map_pool; curvature_map_pool(:);];
            BA_output.filament_alignment_map_pool = [BA_output.filament_alignment_map_pool; filament_alignment_map_pool(:);];
            BA_output.filament_density_pool = [BA_output.filament_density_pool; filament_density_pool(:);];
            BA_output.prot_or_retr_fila_pool = [BA_output.prot_or_retr_fila_pool; prot_or_retr_fila_pool(:);];
            BA_output.int_map_pool = [BA_output.int_map_pool; int_map_pool(:);];
            BA_output.st_map_pool = [BA_output.st_map_pool; st_map_pool(:);];
            BA_output.scale_map_pool = [BA_output.scale_map_pool; scale_map_pool(:);];
            
            BA_output.curvature_map_full_pool = [BA_output.curvature_map_full_pool; curvature_map_full_pool(:);];
            BA_output.filament_alignment_map_full_pool = [BA_output.filament_alignment_map_full_pool; filament_alignment_map_full_pool(:);];
            BA_output.filament_density_full_pool = [BA_output.filament_density_full_pool; filament_density_full_pool(:);];
            BA_output.prot_or_retr_fila_full_pool = [BA_output.prot_or_retr_fila_full_pool; prot_or_retr_fila_full_pool(:);];
            BA_output.int_map_full_pool = [BA_output.int_map_full_pool; int_map_full_pool(:);];
            BA_output.st_map_full_pool = [BA_output.st_map_full_pool; st_map_full_pool(:);];
            
            % what is happening in the previous frame before protusion or
            % retaction, same data, different label as
            % this: have protruded or retracted
            % last: will pretrude or retract to the next time point, so PRE-
            
            last_current_seg = current_seg_cell{1,iCompleteFrame-1};
            
            last_curvature_map = curvature_map_cell{iCompleteFrame};
            last_filament_alignment_map = filament_alignment_map_cell{iCompleteFrame};
            last_filament_density = filament_density_cell{iCompleteFrame};
            last_frame_prot_retract_map_based_on_branch = this_frame_prot_retract_map_based_on_branch_cell{iCompleteFrame};
            
            last_curvature_map_pool = last_curvature_map(last_current_seg>0);
            last_filament_alignment_map_pool = last_filament_alignment_map(last_current_seg>0);
            last_filament_density_pool = last_filament_density(last_current_seg>0);
            last_prot_or_retr_fila_pool = last_frame_prot_retract_map_based_on_branch(last_current_seg>0);
            last_scale_map_pool = last_scaleMap(last_current_seg>0);
            
            last_curvature_map(isnan(last_curvature_map))=0;
            last_filament_alignment_map(isnan(last_filament_alignment_map))=0;
            
            last_curvature_map_full = imfilter(last_curvature_map,ones(radius, radius)/(radius*radius), 'same','replicate')./last_filament_density;
            last_filament_alignment_map_full = imfilter(last_filament_alignment_map,ones(radius, radius)/(radius*radius), 'same','replicate')./last_filament_density;
            
            
            last_curvature_map_full_pool = last_curvature_map_full(last_filament_density>0);
            last_filament_alignment_map_full_pool = last_filament_alignment_map_full(last_filament_density>0);
            last_filament_density_full_pool = last_filament_density(last_filament_density>0);
            last_prot_or_retr_fila_full_pool = last_frame_prot_retract_map_based_on_branch(last_filament_density>0);
            
            
            last_int_map_pool = last_VIF_image(last_current_seg>0);
            last_st_map_pool = this_nms(last_current_seg>0);
            
            last_int_map_full_pool = last_VIF_image(last_filament_density>0);
            last_st_map_full_pool = last_MAX_st_res(last_filament_density>0);
            
            
            BA_output.last_curvature_map_pool = [BA_output.last_curvature_map_pool; last_curvature_map_pool(:);];
            BA_output.last_filament_alignment_map_pool = [BA_output.last_filament_alignment_map_pool; last_filament_alignment_map_pool(:);];
            BA_output.last_filament_density_pool = [BA_output.last_filament_density_pool; last_filament_density_pool(:);];
            BA_output.last_prot_or_retr_fila_pool = [BA_output.last_prot_or_retr_fila_pool; last_prot_or_retr_fila_pool(:);];
            BA_output.last_int_map_pool = [BA_output.last_int_map_pool; last_int_map_pool(:);];
            BA_output.last_st_map_pool = [BA_output.last_st_map_pool; last_st_map_pool(:);];
            BA_output.last_scale_map_pool =[BA_output.last_scale_map_pool; last_scale_map_pool(:);];
            
            BA_output.last_curvature_map_full_pool = [BA_output.last_curvature_map_full_pool; last_curvature_map_full_pool(:);];
            BA_output.last_filament_alignment_map_full_pool = [BA_output.last_filament_alignment_map_full_pool; last_filament_alignment_map_full_pool(:);];
            BA_output.last_filament_density_full_pool = [BA_output.last_filament_density_full_pool; last_filament_density_full_pool(:);];
            BA_output.last_prot_or_retr_fila_full_pool = [BA_output.last_prot_or_retr_fila_full_pool; last_prot_or_retr_fila_full_pool(:);];
            BA_output.last_int_map_full_pool = [BA_output.last_int_map_full_pool; last_int_map_full_pool(:);];
            BA_output.last_st_map_full_pool = [BA_output.last_st_map_full_pool; last_st_map_full_pool(:);];
            
      
            %%
            BA_output.vim_intensity_cell_frame_avg(iCompleteFrame) = mean(current_VIF_image(this_cell_mask>0));
            BA_output.vif_densiity_cell_frame_avg(iCompleteFrame) = mean(filament_density(this_cell_mask>0));
            BA_output.vim_intensity_filaall_frame_avg(iCompleteFrame) = mean(current_VIF_image(current_seg>0));
            BA_output.vif_densiity_filaall_frame_avg(iCompleteFrame) = mean(filament_density(current_seg>0));
            
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
        
        
        
        h15 = figure(15);
        subplot(3,2,1);
        hist(BA_output.curvature_map_pool(BA_output.prot_or_retr_fila_pool==1),30);
        subplot(3,2,2);
        hist(BA_output.curvature_map_pool(BA_output.prot_or_retr_fila_pool==-1),30);
        subplot(3,2,3);
        hist(BA_output.filament_alignment_map_pool(BA_output.prot_or_retr_fila_pool==1),30);
        subplot(3,2,4);
        hist(BA_output.filament_alignment_map_pool(BA_output.prot_or_retr_fila_pool==-1),30);
        subplot(3,2,5);
        hist(BA_output.filament_density_pool(BA_output.prot_or_retr_fila_pool==1),30);
        subplot(3,2,6);
        hist(BA_output.filament_density_pool(BA_output.prot_or_retr_fila_pool==-1),30);
        saveas(h15,[outputPath,filesep,'this_filament_prot_retract_center.tif']);
        saveas(h15,[outputPath,filesep,'this_filament_prot_retract_center.fig']);
        
        
        h16 = figure(16);
        subplot(3,2,1);
        hist(BA_output.curvature_map_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('Curvature full prot-ed');
        subplot(3,2,2);
        hist(BA_output.curvature_map_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('Curvature full retract-ed');
        subplot(3,2,3);
        hist(BA_output.filament_alignment_map_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('Filament Alignment full prot-ed');
        subplot(3,2,4);
        hist(BA_output.filament_alignment_map_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('Filament Alignment full retract-ed');
        subplot(3,2,5);
        hist(BA_output.filament_density_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('Filament density full prot-ed');
        subplot(3,2,6);
        hist(BA_output.filament_density_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('Filament density full retract-ed');
        saveas(h16,[outputPath,filesep,'this_filament_prot_retract_full.tif']);
        saveas(h16,[outputPath,filesep,'this_filament_prot_retract_full.fig']);
        
        h17 = figure(17);
        subplot(3,2,1);
        hist(BA_output.last_curvature_map_full_pool(BA_output.last_prot_or_retr_fila_full_pool==1),30);
        title('Curvature full pre-prot');
        subplot(3,2,2);
        hist(BA_output.last_curvature_map_full_pool(BA_output.last_prot_or_retr_fila_full_pool==-1),30);
        title('Curvature full pre-retract');
        subplot(3,2,3);
        hist(BA_output.last_filament_alignment_map_full_pool(BA_output.last_prot_or_retr_fila_full_pool==1),30);
        title('Filament Alignment full pre-prot');
        subplot(3,2,4);
        hist(BA_output.last_filament_alignment_map_full_pool(BA_output.last_prot_or_retr_fila_full_pool==-1),30);
        title('Filament Alignment full pre-retract');
        subplot(3,2,5);
        hist(BA_output.last_filament_density_full_pool(BA_output.last_prot_or_retr_fila_full_pool==1),30);
        title('Filament density full pre-prot');
        subplot(3,2,6);
        hist(BA_output.last_filament_density_full_pool(BA_output.last_prot_or_retr_fila_full_pool==-1),30);
        title('Filament density full pre-retract');
        
        saveas(h17,[outputPath,filesep,'last_filament_prot_retract_full.tif']);
        saveas(h17,[outputPath,filesep,'last_filament_prot_retract_full.fig']);
        
        h18 = figure(18);
        subplot(3,2,1);
        hist(double(BA_output.int_map_full_pool(BA_output.prot_or_retr_fila_full_pool==1)),30);
        title('Int full prot-ed');
        subplot(3,2,2);
        hist(double(BA_output.int_map_full_pool(BA_output.prot_or_retr_fila_full_pool==-1)),30);
        title('Int full retract-ed');
        subplot(3,2,3);
        hist(double(BA_output.last_int_map_full_pool(BA_output.last_prot_or_retr_fila_full_pool==1)),30);
        title('Int full pre-prot');
        subplot(3,2,4);
        hist(double(BA_output.last_int_map_full_pool(BA_output.last_prot_or_retr_fila_full_pool==-1)),30);
        title('Int full pre-retract');
        subplot(3,2,5);
        hist(BA_output.st_map_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('ST full prot-ed');
        subplot(3,2,6);
        hist(BA_output.st_map_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('ST full retract-ed');
        saveas(h18,[outputPath,filesep,'this_filament_intst_prot_retract_full.tif']);
        saveas(h18,[outputPath,filesep,'this_filament_intst_prot_retract_full.fig']);
        
        
        
        
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
