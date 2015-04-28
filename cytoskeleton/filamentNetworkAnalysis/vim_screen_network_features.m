function output_feature = vim_screen_network_features(labelMaskNucleus,VIF_current_seg,VIF_ROI_model,current_img, nms,feature_flag,T_dis_perp)
% function to calculate network features specifically for vim screen where
% the nucleus segmentation is present
% input:

% VIF_current_seg    the filament segmented
% current_img        the intensity(usually is image-flattened results, but
%                       in the new version, this is no per movie
%                       normalization)
% MAX_st_res         the  steerable filtering results
% nms                the nms version of ST
% T_dis_perp         the division of perphery and center, default 40 pixels

% default
if(nargin<7)
    T_dis_perp = 40;
end

[LabelMask,NoRegion] = bwlabel(labelMaskNucleus);

NewLabel = randperm(NoRegion);
NewLabelMask = LabelMask;
for iR = 1 :NoRegion
    NewLabelMask(LabelMask==iR) = NewLabel(iR);
end

Mask = labelMaskNucleus;

[DistMask,IDX] = bwdist(Mask);
RegionMask = NewLabelMask;
RegionMask(:) = NewLabelMask(IDX);

nucleus_center_x = [];
nucleus_center_y = [];
for  iR = 1 : NoRegion
    [ind_y, ind_x ] = find(NewLabelMask==iR);
    nucleus_center_x(iR)=mean(ind_x);
    nucleus_center_y(iR)=mean(ind_y);
end

%%
if(feature_flag(19)>0 || feature_flag(20)>0)
    
    profileCell = cell(1,NoRegion);
    
    display_regions_sections = NewLabelMask;
    
    for  iR = 1 : NoRegion
        region_thisCell = RegionMask==iR;
        nucleus_thisCell = NewLabelMask==iR;
        distMap_thisCell = DistMask;
        
        % to include the boundary of the nucleus
        shrinked_nucleus_thisCell = imerode(nucleus_thisCell,[ 0 1 0; 1 1 1; 0 1 0]);
        
        inverted_Distmap = (bwdist(1-shrinked_nucleus_thisCell));
        % get the maximum points
        max_invertdist = inverted_Distmap==max(max(inverted_Distmap));
        
        % in case there are more than one of these max centers
        max_invertdist = keep_largest_area(max_invertdist);
        
        % in each maximum, there could be more than one point, to make the
        % code easier for its own flow
        C_xy = regionprops(max_invertdist,'Centroid');
        center_x = C_xy.Centroid(1);
        center_y = C_xy.Centroid(2);
        
        region_without_shrinked_nucleus_thisCell = region_thisCell - shrinked_nucleus_thisCell;
        
        distMap_thisCell = distMap_thisCell(region_without_shrinked_nucleus_thisCell>0);
        VIF_current_seg_thisCell = VIF_current_seg(region_without_shrinked_nucleus_thisCell>0);
        current_img_thisCell = current_img(region_without_shrinked_nucleus_thisCell>0);
        stnms_thisCell =  nms(region_without_shrinked_nucleus_thisCell>0);
        
        % Checking the distance map,that is the distance from the nucleus
        % to get the profile of the vim/mt
        % intensity/filamentdensity  as the distance increase from 0( close to
        % nucleus) to max(a quarter of the image size(1))
        
        [indy,indx] = find(region_without_shrinked_nucleus_thisCell>0);
        
        angle_to_nucleus_center = atan2((indy-center_y),(indx-center_x));
        
        profileIntSum = nan(6,max(distMap_thisCell)+1);
        profileFilaSum = nan(6,max(distMap_thisCell)+1);        
        profileStnmsSum = nan(6,max(distMap_thisCell)+1);         
        profilePixels = nan(6,max(distMap_thisCell)+1);        
%         profileIntOverST = nan(6,max(distMap_thisCell)+1);       

        for iDistance = 1: max(distMap_thisCell)+1;
            DistanceT = iDistance-1;
            for iAngle =  1 : 6
                Angle_bottom = -pi + (2*pi/6) *(iAngle-1);
                Angle_top =  -pi + (2*pi/6) *(iAngle);
                
                ind_this_section  = find(round(distMap_thisCell)==DistanceT & ...
                    angle_to_nucleus_center>=Angle_bottom & angle_to_nucleus_center<=Angle_top);
                
                display_regions_sections(...
                    sub2ind(size(NewLabelMask),indy(ind_this_section),indx(ind_this_section)))...
                    =iR-0.5+iAngle/7;
                
                profileIntSum(iAngle, iDistance) = ...
                    sum(current_img_thisCell(ind_this_section));
                profileFilaSum(iAngle, iDistance) = ...
                    sum(VIF_current_seg_thisCell(ind_this_section));
                profileStnmsSum(iAngle, iDistance) =  ...
                    sum(stnms_thisCell(ind_this_section));
                profilePixels(iAngle, iDistance) = numel(ind_this_section);
                
                
            end
            
        end
        
%         DistancePeriCenter = nan(1,6);
%         
%         for iAngle =  1 : 6
%             
%             profileIntAve_array = profileIntSum(iAngle, :)./profilePixels(iAngle, :);
%             profileSTAve_array = profileStnmsSum(iAngle, :)./profilePixels(iAngle, :);
%             
%             profileIntAve_array_smooth = imfilter(profileIntAve_array, fspecial('gaussian', 11, 4), 'replicate','same');
%             profileSTAve_array_smooth = imfilter(profileSTAve_array, fspecial('gaussian', 11, 4), 'replicate','same');
%             
%             STOverInt_array = profileSTAve_array_smooth./profileIntAve_array_smooth;            
%             STOverInt_array_smooth = imfilter(STOverInt_array, fspecial('gaussian', 11, 4), 'replicate','same');
%             STOverInt_array_grad = STOverInt_array_smooth(2:end)-STOverInt_array_smooth(1:end-1);
%             
%             ind_d = find(STOverInt_array_grad == max(STOverInt_array_grad));
%             DistancePeriCenter(iAngle) = ind_d;
%         end
        
        profileCell{1,iR}.profileIntSum = profileIntSum;
        profileCell{1,iR}.profileFilaSum = profileFilaSum;
        profileCell{1,iR}.profileStnmsSum = profileStnmsSum;
        profileCell{1,iR}.profilePixels = profilePixels;
        profileCell{1,iR}.profileIntMean = profileIntSum./profilePixels;
        profileCell{1,iR}.profileFilaMean = profileFilaSum./profilePixels;
        profileCell{1,iR}.profileStnmsMean = profileStnmsSum./profilePixels;
        
        try
            profileCell{1,iR}.profileIntSumPerpCenterRatio = ...
                sum(profileIntSum(1:6,T_dis_perp:end),2)./sum(profileIntSum(1:6,1:T_dis_perp+1),2);
        catch
            profileCell{1,iR}.profileIntSumPerpCenterRatio = nan;
        end
        
        try
            profileCell{1,iR}.profileIntMeanPerpCenterRatio = ...
                (sum(profileIntSum(1:6,T_dis_perp+1:end),2)./sum(profilePixels(1:6,T_dis_perp+1:end),2))...
                ./(sum(profileIntSum(1:6,1:T_dis_perp),2)./sum(profilePixels(1:6,1:T_dis_perp),2));
        catch
            profileCell{1,iR}.profileIntMeanPerpCenterRatio = nan;
        end
        
        try            
            profileCell{1,iR}.profileFilaSumPerpCenterRatio = ...
                sum(profileFilaSum(1:6,T_dis_perp+1:end),2)...
                ./sum(profileFilaSum(1:6,1:T_dis_perp),2);
        catch
            profileCell{1,iR}.profileFilaSumPerpCenterRatio = nan;
        end
        
        try            
            profileCell{1,iR}.profileFilaMeanPerpCenterRatio = ...
                (sum(profileFilaSum(1:6,T_dis_perp+1:end),2)./sum(profilePixels(1:6,T_dis_perp+1:end),2))...
                ./(sum(profileFilaSum(1:6,1:T_dis_perp),2)./sum(profilePixels(1:6,1:T_dis_perp),2));
        catch
            profileCell{1,iR}.profileFilaMeanPerpCenterRatio = nan;
        end
        
        try            
            profileCell{1,iR}.profileStnmsSumPerpCenterRatio = ...
                sum(profileStnmsSum(1:6,T_dis_perp+1:end),2)...
                ./sum(profileStnmsSum(1:6,1:T_dis_perp),2);
        catch
            profileCell{1,iR}.profileStnmsSumPerpCenterRatio = nan;
        end
        
        try            
            profileCell{1,iR}.profileStnmsMeanPerpCenterRatio = ...
                (sum(profileStnmsSum(1:6,T_dis_perp+1:end),2)./sum(profilePixels(1:6,T_dis_perp+1:end),2))...
                ./(sum(profileStnmsSum(1:6,1:T_dis_perp),2)./sum(profilePixels(1:6,1:T_dis_perp),2));
        catch
            profileCell{1,iR}.profileStnmsMeanPerpCenterRatio = nan;
        end
    end
    
    intsum_pool = [];
    intmean_pool = [];
    filasum_pool = [];
    filamean_pool = [];
    nmssum_pool = [];
    nmsmean_pool = [];
    
    for  iR = 1 : NoRegion
        nmssum_pool = [nmssum_pool profileCell{1,iR}.profileStnmsSumPerpCenterRatio];
        nmsmean_pool = [nmsmean_pool profileCell{1,iR}.profileStnmsMeanPerpCenterRatio];
        intsum_pool = [intsum_pool profileCell{1,iR}.profileIntSumPerpCenterRatio];
        intmean_pool = [intmean_pool profileCell{1,iR}.profileIntMeanPerpCenterRatio];
        filasum_pool = [filasum_pool profileCell{1,iR}.profileFilaSumPerpCenterRatio];
        filamean_pool = [filamean_pool profileCell{1,iR}.profileFilaMeanPerpCenterRatio];
    end
    
    profileAllCell.MedStnmsSumPerpCenterRatio  = nanmedian(nmssum_pool(:));
    profileAllCell.MedStnmsMeanPerpCenterRatio = nanmedian(nmsmean_pool(:));
    profileAllCell.MedIntSumPerpCenterRatio    = nanmedian(intsum_pool(:));
    profileAllCell.MedIntMeanPerpCenterRatio   = nanmedian(intmean_pool(:));
    profileAllCell.MedFilaSumPerpCenterRatio   = nanmedian(filasum_pool(:));
    profileAllCell.MedFilaMeanPerpCenterRatio  = nanmedian(filamean_pool(:));
    
    profileAllCell.MeanStnmsSumPerpCenterRatio  = nanmean(nmssum_pool(:));
    profileAllCell.MeanStnmsMeanPerpCenterRatio = nanmean(nmsmean_pool(:));
    profileAllCell.MeanIntSumPerpCenterRatio    = nanmean(intsum_pool(:));
    profileAllCell.MeanIntMeanPerpCenterRatio   = nanmean(intmean_pool(:));
    profileAllCell.MeanFilaSumPerpCenterRatio   = nanmean(filasum_pool(:));
    profileAllCell.MeanFilaMeanPerpCenterRatio  = nanmean(filamean_pool(:));
    
    profileAllCell.nmssum_pool = nmssum_pool;
    profileAllCell.nmsmean_pool = nmsmean_pool;
    profileAllCell.intsum_pool = intmean_pool;
    profileAllCell.intmean_pool = intmean_pool;
    profileAllCell.filasum_pool = filasum_pool;
    profileAllCell.filamean_pool = filamean_pool;
    
     output_feature.profileCell = profileCell;
    output_feature.profileAllCell = profileAllCell;
else    
     output_feature.profileCell = [];
    output_feature.profileAllCell = [];    
end
%% Define centripetal-ness
if(feature_flag(21)>0 ||feature_flag(22)>0)
    angle_array=[];
    angle_pixel_array=[];
%     figure(20);imagesc(NewLabelMask);
%     CCCCC = colormap;
%     CCCCC(1,:)=[1 1 1];
%     colormap(CCCCC);
    
    for iF = 1 : length(VIF_ROI_model)
        
        x = VIF_ROI_model{iF}(:,1);
        y = VIF_ROI_model{iF}(:,2);
        
        x = x(:);
        y = y(:);
        
        smooth_x = imfilter(x, [ 1 2 5 2 1]'/11,'same','replicate');
        smooth_y = imfilter(y, [ 1 2 5 2 1]'/11,'same','replicate');
        
        orientation_pixel_array = atan2(smooth_y(2:end)-smooth_y(1:end-1),...
            smooth_x(2:end)-smooth_x(1:end-1));
        
        F_x = x(round(numel(x)/2));
        F_y = y(round(numel(x)/2));
        
        orientation_filament = orientation_pixel_array(round(numel(x)/2));
        
        region_ind_array =  RegionMask(sub2ind(size(NewLabelMask),y, x));
        region_ind = mode(region_ind_array);
        
        N_x = nucleus_center_x(region_ind);
        N_y = nucleus_center_y(region_ind);
        
        orientation_connecting = atan2(F_y - N_y, F_x - N_x);
        angle_array(iF) = orientation_connecting - orientation_filament;
        
%         hold on;
%         this_color_array = CCCCC(round(size(CCCCC,1)/NoRegion*region_ind),:);
%         hold on; plot(smooth_x,smooth_y,'-','color',this_color_array/1.1);
%         
%         line_this_color_array = this_color_array*1.2;
%         line_this_color_array(line_this_color_array>1)=1;
%         hold on; plot([N_x F_x],[N_y F_y],'-','color',line_this_color_array);
       
        % Liya, the pixel based should be more reasonable
        
        orientation_pixel_connecting = atan2(y(1:end-1) - N_y, x(1:end-1) - N_x);
       
        angle_pixel_array = [angle_pixel_array;...
            (orientation_pixel_array - orientation_pixel_connecting)];
        
            end
    
    % mod to 0~pi
    angle_array(angle_array>pi) = angle_array(angle_array>pi) -pi;
    angle_array(angle_array<0) = angle_array(angle_array<0) +pi;
    angle_array(angle_array>pi) = angle_array(angle_array>pi) -pi;
    angle_array(angle_array<0) = angle_array(angle_array<0) +pi;
    % now everything 0~pi, then convert to 0~pi/2
    angle_array(angle_array>pi/2) = pi - angle_array(angle_array>pi/2) ;
    
%     HH = histc(angle_array,0:pi/24:pi/2);
%     % combine the ==pi/2 into the last valid bin
%     HH(end-1) = HH(end-1) +HH(end);
%     HH = HH(1:end-1);
%     
%     Centripetal_distribution = HH;
    output_feature.Centripetal_fila = angle_array;
    
    % pixel based
    
    % mod to 0~pi
    angle_pixel_array(angle_pixel_array>pi) = angle_pixel_array(angle_pixel_array>pi) -pi;
    angle_pixel_array(angle_pixel_array<0) = angle_pixel_array(angle_pixel_array<0) +pi;
    angle_pixel_array(angle_pixel_array>pi) = angle_pixel_array(angle_pixel_array>pi) -pi;
    angle_pixel_array(angle_pixel_array<0) = angle_pixel_array(angle_pixel_array<0) +pi;
    % now everything 0~pi, then convert to 0~pi/2
    angle_pixel_array(angle_pixel_array>pi/2) = pi - angle_pixel_array(angle_pixel_array>pi/2) ;
    
%     HH = histc(angle_pixel_array,0:pi/24:pi/2);
%     % combine the ==pi/2 into the last valid bin
%     HH(end-1) = HH(end-1) +HH(end);
%     HH = HH(1:end-1);
%     
%     Centripetal_pixel_distribution = HH;
    output_feature.Centripetal_pixel = angle_pixel_array;
else
    output_feature.Centripetal_fila = [];
    output_feature.Centripetal_pixel = [];
end
