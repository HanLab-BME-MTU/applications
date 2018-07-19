 
%         
%   [MT_digital_model,MT_orientation_model,MT_XX,MT_YY,MT_OO] ...
%         = filament_model_to_digital_with_orientation(MT_current_model);
%   [VIF_digital_model,VIF_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
%         = filament_model_to_digital_with_orientation(VIF_current_model);

    
    
    
        Hue = (VIF_OO(:)+pi/2)/(pi)-0.2;
        Hue(find(Hue>=1)) = Hue(find(Hue>=1)) -1;
        Hue(find(Hue<0)) = Hue(find(Hue<0)) +1;
        
        Sat = Hue*0+1;
        Value = Hue*0+1;
       
        RGB_seg_orient_heat_array = hsv2rgb([Hue Sat Value]);
        
        R_seg_orient_heat_map = nan(img_size);
        G_seg_orient_heat_map = nan(img_size);
        B_seg_orient_heat_map = nan(img_size);
        
        R_seg_orient_heat_map(sub2ind(img_size,VIF_YY,VIF_XX))=RGB_seg_orient_heat_array(:,1);
        G_seg_orient_heat_map(sub2ind(img_size,VIF_YY,VIF_XX))=RGB_seg_orient_heat_array(:,2);
        B_seg_orient_heat_map(sub2ind(img_size,VIF_YY,VIF_XX))=RGB_seg_orient_heat_array(:,3);
        
        RGB_seg_orient_heat_map = nan(img_size(1),img_size(2),3);
        
        RGB_seg_orient_heat_map(:,:,1) = R_seg_orient_heat_map;
        RGB_seg_orient_heat_map(:,:,2) = G_seg_orient_heat_map;
        RGB_seg_orient_heat_map(:,:,3) = B_seg_orient_heat_map;
        
        figure;imagesc(RGB_seg_orient_heat_map);
        hold on;
        plot(X1(crossing_flag_1>0),Y1(crossing_flag_1>0),'o','MarkerSize',10,...
            'MarkerFaceColor','g','MarkerEdgeColor','k');
     plot(X1(crossing_flag_1>0),Y1(crossing_flag_1>0),'rx','MarkerSize',12,...
            'MarkerFaceColor','r');
    
        
        %%
        
%         VIF_current_seg
%         VIF_orientation
        
        [Y_seg,X_seg] = find(VIF_current_seg>0);
        O_seg = VIF_orientation(sub2ind(img_size,Y_seg,X_seg))
        
        
        Hue = (-O_seg(:)+pi/2)/(pi)-0.2;
        Hue(find(Hue>=1)) = Hue(find(Hue>=1)) -1;
        Hue(find(Hue<0)) = Hue(find(Hue<0)) +1;
        
        Sat = Hue*0+1;
        Value = Hue*0+1;
       
        RGB_seg_orient_heat_array = hsv2rgb([Hue Sat Value]);
        
        R_seg_orient_heat_map = nan(img_size);
        G_seg_orient_heat_map = nan(img_size);
        B_seg_orient_heat_map = nan(img_size);
        
        R_seg_orient_heat_map(sub2ind(img_size,Y_seg,X_seg))=RGB_seg_orient_heat_array(:,1);
        G_seg_orient_heat_map(sub2ind(img_size,Y_seg,X_seg))=RGB_seg_orient_heat_array(:,2);
        B_seg_orient_heat_map(sub2ind(img_size,Y_seg,X_seg))=RGB_seg_orient_heat_array(:,3);
        
        RGB_seg_orient_heat_map_seg = nan(img_size(1),img_size(2),3);
        
        RGB_seg_orient_heat_map_seg(:,:,1) = R_seg_orient_heat_map;
        RGB_seg_orient_heat_map_seg(:,:,2) = G_seg_orient_heat_map;
        RGB_seg_orient_heat_map_seg(:,:,3) = B_seg_orient_heat_map;
        
        figure;imagesc(RGB_seg_orient_heat_map_seg);
        
        