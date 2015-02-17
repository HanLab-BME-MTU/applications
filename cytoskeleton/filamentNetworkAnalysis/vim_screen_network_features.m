function output_feature = vim_screen_network_features(labelMaskNucleus,VIF_current_seg,current_img, nms)
% function to calculate network features specifically for vim screen where
% the nucleus segmentation is present
% input:

% VIF_current_seg    the filament segmented
% current_img        the intensity(usually is image-flattened results, but
%                       in the new version, this is no per movie
%                       normalization)
% MAX_st_res         the  steerable filtering results
% nms                the nms version of ST
% 

[LabelMask,NoRegion] = bwlabel(labelMaskNucleus);

NewLabel = randperm(NoRegion);
NewLabelMask = LabelMask;
for iR = 1 :NoRegion
    NewLabelMask(LabelMask==iR) = NewLabel(iR);
end

[DistMask,IDX] = bwdist(Mask);
RegionMask = NewLabelMask;
RegionMask(:) = NewLabelMask(IDX);

profileCell = cell(1,NoRegion);

for  iR = 1 : NoRegion
    region_thisCell = RegionMask==iR;
    nucleus_thisCell = NewLabelMask==iR;
    
    % to include the boundary of the nucleus
    shrinked_nucleus_thisCell = imerode(nucleus_thisCell);
    
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
    
    for iDistance = 1: max(distMap_thisCell)+1;
        DistanceT = iDistance-1;
        for iAngle =  1 : 6
            Angle_bottom = -pi + (2*pi/6) *(iAngle-1);
            Angle_top =  -pi + (2*pi/6) *(iAngle);
            
            ind_this_section  = find(round(distMap_thisCell)==DistanceT & ...
                angle_to_nucleus_center>=Angle_bottom & angle_to_nucleus_center<=Angle_top);
            
            profileIntSum(iAngle, iDistance) = ...
                sum(current_img_thisCell(ind_this_section));
            profileFilaSum(iAngle, iDistance) = ...
                sum(VIF_current_seg_thisCell(ind_this_section));
            profileStnmsSum(iAngle, iDistance) =  ...
                sum(stnms_thisCell(ind_this_section));
            profilePixels(iAngle, iDistance) = numel(ind_this_section);
        end
    end
    
    profileCell{1,iR}.profileIntSum = profileIntSum;
    profileCell{1,iR}.profileFilaSum = profileFilaSum;
    profileCell{1,iR}.profileStnmsSum = profileStnmsSum;
    profileCell{1,iR}.profilePixels = profilePixels;
    
    %% Define centripetal-ness
    VIF_ROI_model
    
    for iF = 1 : length(VIF_ROI_model)
        try
            x = VIF_ROI_model{iF}(:,1);
            y = VIF_ROI_model{iF}(:,2);
    
    
end

output_feature.profileCell = profileCell;    