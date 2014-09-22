function [filament_similiarity_1,filament_similiarity_2,filament_similarity_scoremap] ...
    = filament_similarity(VIF_current_model,MT_current_model,img_size, ...
    radius,outdir,iFrame,save_everything_flag,...
    distance_map_1_2, distance_map_2_1, angle_map_1_2, angle_map_2_1)
% function for calculation the similarity of two networks
% on the base of filament by filament
% Liya Ding 05.2014.

% don't save every figure generated, unless debugging
% save_everything_flag = 0

[MT_digital_model,MT_orientation_model,MT_XX,MT_YY,MT_OO] ...
    = filament_model_to_digital_with_orientation(MT_current_model);
[VIF_digital_model,VIF_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
    = filament_model_to_digital_with_orientation(VIF_current_model);

filament_similiarity_1=[];
filament_similiarity_2=[];

filament_similiarity_1.distance_map_1_2_filament = distance_map_1_2;
filament_similiarity_2.distance_map_2_1_filament = distance_map_2_1;
filament_similiarity_1.angle_map_1_2_filament = angle_map_1_2;
filament_similiarity_2.angle_map_2_1_filament = angle_map_2_1;

nF_1 = length(MT_digital_model);
for iF = 1 : nF_1
    % find the value for all pixels on this filament
    DD = distance_map_1_2(sub2ind(img_size,MT_digital_model{iF}(:,2),MT_digital_model{iF}(:,1)));
    AA = angle_map_1_2(sub2ind(img_size,MT_digital_model{iF}(:,2),MT_digital_model{iF}(:,1)));
    
    % take the nan mean as the filament feature
    filament_similiarity_1.distance(iF) = nanmean(DD);
    filament_similiarity_1.angle(iF) = nanmean(AA);
    
    % replace the pixel value with the filament mean
    filament_similiarity_1.distance_map_1_2_filament(sub2ind(img_size,MT_digital_model{iF}(:,2),MT_digital_model{iF}(:,1)))=filament_similiarity_1.distance(iF);
    filament_similiarity_1.angle_map_1_2_filament(sub2ind(img_size,MT_digital_model{iF}(:,2),MT_digital_model{iF}(:,1)))=filament_similiarity_1.angle(iF);
    
end

nF_2 = length(VIF_digital_model);
for iF = 1 : nF_2
    DD = distance_map_2_1(sub2ind(img_size,VIF_digital_model{iF}(:,2),VIF_digital_model{iF}(:,1)));
    AA = angle_map_2_1(sub2ind(img_size,VIF_digital_model{iF}(:,2),VIF_digital_model{iF}(:,1)));
    
    filament_similiarity_2.distance(iF) = nanmean(DD);
    filament_similiarity_2.angle(iF) = nanmean(AA);
    
    filament_similiarity_2.distance_map_2_1_filament(sub2ind(img_size,VIF_digital_model{iF}(:,2),VIF_digital_model{iF}(:,1)))=filament_similiarity_2.distance(iF);
    filament_similiarity_2.angle_map_2_1_filament(sub2ind(img_size,VIF_digital_model{iF}(:,2),VIF_digital_model{iF}(:,1)))=filament_similiarity_2.angle(iF);
    
end

%% build a new scoremap based on updated difference

distance_map_1_2_pad = nan(img_size+2*radius);
distance_map_2_1_pad = nan(img_size+2*radius);
angle_map_1_2_pad = nan(img_size+2*radius);
angle_map_2_1_pad = nan(img_size+2*radius);


distance_map_1_2_pad(radius+1:end-radius,radius+1:end-radius) = filament_similiarity_1.distance_map_1_2_filament;
distance_map_2_1_pad(radius+1:end-radius,radius+1:end-radius) = filament_similiarity_2.distance_map_2_1_filament;
angle_map_1_2_pad(radius+1:end-radius,radius+1:end-radius) = filament_similiarity_1.angle_map_1_2_filament;
angle_map_2_1_pad(radius+1:end-radius,radius+1:end-radius) = filament_similiarity_2.angle_map_2_1_filament;

% size of the disk matrix is 2*radius+1
[cy,cx] = find(fspecial('disk',radius)>0);

% size of the kernal is also 2*radius+1
Weight_mask = fspecial('gaussian',2*radius+1,(radius*1.3)/(2)/2);

% make the matrix into index array
Weight_mask = Weight_mask(sub2ind([2*radius+1,2*radius+1,],cy,cx));

filament_score_maps_distance_1_2 = nan(img_size);
filament_score_maps_distance_2_1 = nan(img_size);

filament_score_maps_angle_1_2 = nan(img_size);
filament_score_maps_angle_2_1 = nan(img_size);


VIF_current_seg = filament_model_to_seg_bwim(VIF_current_model,img_size,[]);
MT_current_seg = filament_model_to_seg_bwim(MT_current_model,img_size,[]);
whole_ROI = imdilate(VIF_current_seg,ones(radius,radius)) + imdilate(MT_current_seg,ones(radius,radius))>0;
% find the points of interest
[Y,X] = find(whole_ROI>0);

% for all these points, get local support
for j = 1 : length(Y)
%     j
    x = X(j);
    y = Y(j);
        
    dis_1 = distance_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));
    dis_2 = distance_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));

    ang_1 = abs(angle_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1)));
    ang_2 = abs(angle_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1)));
    
    if(~isempty(dis_1))
        filament_score_maps_distance_1_2(y,x) = sum(dis_1(~isnan(dis_1)).*Weight_mask(~isnan(dis_1)))/sum(Weight_mask(~isnan(dis_1)));
        filament_score_maps_angle_1_2(y,x) = sum(ang_1(~isnan(ang_1)).*Weight_mask(~isnan(ang_1)))/sum(Weight_mask(~isnan(ang_1)));
    end
    if(~isempty(dis_2))
        filament_score_maps_distance_2_1(y,x) = sum(dis_2(~isnan(dis_2)).*Weight_mask(~isnan(dis_2)))/sum(Weight_mask(~isnan(dis_2)));       
        filament_score_maps_angle_2_1(y,x) = sum(ang_2(~isnan(ang_2)).*Weight_mask(~isnan(ang_2)))/sum(Weight_mask(~isnan(ang_2)));
    end    
end

% calculation of similarity score.
filament_similarity_scoremap = exp(-(filament_score_maps_distance_2_1+filament_score_maps_distance_1_2).^2/(((radius*1.5)/2*sqrt(2))^2))...
    .*exp(-(abs(filament_score_maps_angle_2_1/2)+abs(filament_score_maps_angle_1_2/2)).^2/(1.5*(pi/3)^2));

% calculation of similarity score only consider 1->2.
filament_similarity_scoremap_1to2 = exp(-(0+filament_score_maps_distance_1_2*2).^2/(((radius*1.5)/2*sqrt(2))^2))...
    .*exp(-(abs(0/2)+abs(filament_score_maps_angle_1_2/2*2)).^2/(1.5*(pi/3)^2));

% calculation of similarity score only consider 2->1.
filament_similarity_scoremap_2to1 = exp(-(filament_score_maps_distance_2_1*2+0).^2/(((radius*1.5)/2*sqrt(2))^2))...
    .*exp(-(abs(filament_score_maps_angle_2_1/2*2)+abs(0/2)).^2/(1.5*(pi/3)^2));

filament_similiarity_1.filament_similarity_scoremap_1to2=filament_similarity_scoremap_1to2;
filament_similiarity_2.filament_similarity_scoremap_2to1=filament_similarity_scoremap_2to1;
filament_similiarity_1.filament_score_maps_distance_1_2=filament_score_maps_distance_1_2;
filament_similiarity_2.filament_score_maps_distance_2_1=filament_score_maps_distance_2_1;
filament_similiarity_1.filament_score_maps_angle_1_2=filament_score_maps_angle_1_2;
filament_similiarity_2.filament_score_maps_angle_2_1=filament_score_maps_angle_2_1;

save([outdir,filesep,'VIFMT_sm_maps_filament_frame_',num2str(iFrame),'.mat'], ...
    'filament_similiarity_1', 'filament_similiarity_2',...
    'filament_similarity_scoremap');

h6=figure(6); imagesc_nan_neg(filament_similarity_scoremap,0);axis equal;axis off;
title(['Filament Similarity Score for frame ',num2str(iFrame)]);

if(save_everything_flag>0)
    saveas(h6,[outdir,filesep,'VIFMT_sm_score_filament_frame_',num2str(iFrame),'.tif']);
    saveas(h6,[outdir,filesep,'VIFMT_sm_score_filament_frame_',num2str(iFrame),'.fig']);
end
