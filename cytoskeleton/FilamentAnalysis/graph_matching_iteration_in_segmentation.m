if retrain_flag==1
    train_mat = [];
    for T_xie_int_grid = T_xie_int_train*(0.8) : (T_xie_int_train*(1.2) - T_xie_int_train*(0.8))/20 : T_xie_int_train*(1.2)
        for T_xie_length_grid = T_xie_length_train*(0.8) : (T_xie_length_train*(1.2) - T_xie_length_train*(0.8))/20 : T_xie_length_train*(1.2)
            
            F_classifer_train = @(i,l) (((T_xie_int_grid + (T_xie_int_grid/T_xie_length_grid)*(-l) -i )));
            train_mat = [train_mat; T_xie_int_grid T_xie_length_grid ...
                (-sum(F_classifer_train(feature_MeanNMS(Matched_ind),...
                feature_Length((Matched_ind))))...
                +sum(F_classifer_train(feature_MeanNMS(UnMatched_ind),...
                feature_Length((UnMatched_ind)))))];
        end
    end
    
    ind = find(train_mat(:,3)==max(train_mat(:,3)));
    
    T_xie_int_train = train_mat(ind(1), 1);
    T_xie_length_train = train_mat(ind(1), 2);
else
    T_xie_int_train = T_xie_int_train*(0.95);
    T_xie_length_train = T_xie_length_train*(0.95);
end

F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
% Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
% Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length)==0);

Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0 & feature_Curvature<CurvatureThreshold & feature_Length>=LengthThreshold);
Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length)==0 | feature_Curvature>=CurvatureThreshold | feature_Length<LengthThreshold);


Good_ind_cell{iIteration+1} = Good_ind;
Bad_ind_cell{iIteration+1} = Bad_ind;

h12 = figure(12);set(h12,'Visible',set_visible);hold off;
plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');hold on;
plot(feature_Length(Matched_ind),feature_MeanNMS(Matched_ind),'g.');
plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.');


title(['Classifier Plane with matched data before round ',num2str(iIteration)]);
if(SaveFigures==1)
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.tif']);
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.fig']);
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.eps']);
end


% plot the output image with these good ones
current_all_seg_bw = zeros(size(labelMask));
for i_E = 1 : length(Good_ind)
    current_good_bw = labelMask==Good_ind(i_E);
    current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
end

imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_begin.tif']);
imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_begin.bmp']);


if(exist('h1','var'))
    close(h1);
end

[current_model,current_matching_bw, model_ind,transformer, final_set_rescue, final_set_connection, connect_count,bw_rgb_three_color,bw_whitergb_three_color]...
    = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt, ...
    Good_ind,ind_long, model_ind,feature_all,labelMask,master_flassier,iIteration,funParams,final_set_rescue, final_set_connection, connect_count,Original_set_Good_ind);

 imwrite(double(bw_whitergb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor_pixel.tif']);
  
if(SaveFigures==1)
    imwrite(double(current_matching_bw/2),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_end.tif']);
    imwrite(double(bw_rgb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor_pixel.tif']);
     
    h1=figure(1);set(h1,'Visible',set_visible);
    saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.tif']);
    saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.fig']);
    print(h1,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.eps']);
    h3=figure(3);set(h3,'Visible',set_visible);
    saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.tif']);
    saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.fig']);
    print(h3,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.eps']);
    
    h11=figure(11);set(h11,'Visible',set_visible);
    saveas(h11,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.tif']);
    saveas(h11,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.fig']);
    print(h11,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.eps']);
    
    h21=figure(21);set(h21,'Visible',set_visible);
    saveas(h21,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.tif']);
    saveas(h21,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.fig']);
    print(h21,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.eps']);
    
end

good_bw = nms_seg_no_brancing.*current_matching_bw;

original_Matched_ind = Matched_ind;

Matched_ind = model_ind(find(~isnan(model_ind)));

Matched_ind_cell{iIteration+1} = Matched_ind;
UnMatched_ind_cell{iIteration+1} = setdiff(ind_long,Matched_ind);
model_ind_cell{iIteration+1} = model_ind;

Good_ind = find(master_flassier(feature_MeanNMS, feature_Length,feature_MeanInt, feature_Curvature)>0);
Bad_ind = find(master_flassier(feature_MeanNMS, feature_Length,feature_MeanInt, feature_Curvature)==0);

h12 = figure(12);set(h12,'Visible',set_visible);hold off;
plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');hold on;
plot(feature_Length(Matched_ind),feature_MeanNMS(Matched_ind),'g.');
plot(feature_Length(setdiff(original_Matched_ind,Good_ind)),feature_MeanNMS(setdiff(original_Matched_ind,Good_ind)),'m*');
plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.');

title(['Classifier Plane with matched data after round ',num2str(iIteration)]);
% saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.tif']);

if(SaveFigures==1)
    print(h12,'-depsc2',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.eps']);
    print(h12,'-dtiff',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.tiff']);
    
    
    for iiii = setdiff(ind_long,UnMatched_ind_cell{iIteration+1})
        text(feature_Length(iiii),feature_MeanNMS(iiii),num2str(iiii));
    end
    
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_',num2str(iFrame),num2str(iIteration),'.fig']);
    
    h13 = figure(13);set(h13,'Visible',set_visible);
    RGB_seg_orient_heat_map = heat_display_filament_from_model(imageInt, current_model);
    imshow(RGB_seg_orient_heat_map);
    title(['Heat display for Segmentation after round ',num2str(iIteration)]);
    imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/GEO/heat_frame_',num2str(iFrame),'_afterround',num2str(iIteration),'.tif']);
    imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/GEO/heat_frame_',num2str(iFrame),'_afterround',num2str(iIteration),'.bmp']);
end
% new_curves = zeros(size(imageNMS));
%
% for ind = (setdiff(Matched_ind,original_Matched_ind))'
% new_curves = imageNMS.*( labelMask==ind);
% title(['ind=',num2str(ind)]);
% % new_curves(1,1)=max(max(imageNMS));
% figure(1);imagescc(new_curves)
% h1 = figure(1);print(h1,'-dtiff',[num2str(ind),'b.tif']);
%
% end

% new_curves = zeros(size(imageNMS));
%
% for ind = (setdiff(Matched_ind,original_Matched_ind))'
% new_curves = new_curves+ imageNMS.*( labelMask==ind);
%
% end
% % new_curves(1,1)=max(max(imageNMS));
% figure(1);imagescc(new_curves)
% h1 = figure(1);print(h1,'-dtiff','brond3_new_all.tif');


