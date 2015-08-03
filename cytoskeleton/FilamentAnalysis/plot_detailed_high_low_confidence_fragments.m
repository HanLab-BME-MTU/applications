h22 =  figure(22);hold off;
plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'.','color',[100 100 255]/255); hold on;
plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');
plot(feature_Length(changed_ind_good_K),feature_MeanNMS(changed_ind_good_K),'c.');
plot(feature_Length(changed_ind_good_I),feature_MeanNMS(changed_ind_good_I),'.','color',[255 128 0   ] / 255);

plot_length = -1:0.1:100;

plot_nms = (T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-plot_length) );

hold on; plot(plot_length, plot_nms,'color',[170 170 170] / 255,'LineWidth',2);
axis([0 60 0 10]);
title('Classifier Plane with matched data round 0');
title('');
set(gca,'FontSize',13)
saveas(h22,['./GEO_frame_',num2str(iFrame),'_round0_trained_plane.tif']);
saveas(h22,['./GEO_frame_',num2str(iFrame),'_round0_trained_plane.fig']);

 
h23 = figure(23);hold off;
imagesc(imageInt_nohole); colormap(gray); axis image; hold on;
axis off;
h26 = figure(26);hold off;
imagesc(imageInt_nohole); colormap(gray); axis image; hold on;
axis off;


% get the mean intensity of the curves
for i_area = setdiff(setdiff(intersect(Good_ind',ind_long'),changed_ind_good_K),changed_ind_good_I)
    [all_y_i, all_x_i] = find(labelMask == i_area);
    % only for those good ones, check the curvature
    
    bw_i = zeros(size(bw_out));
    bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
    end_points_i = bwmorph(bw_i,'endpoints');
    [y_i, x_i]=find(end_points_i);
    
    if isempty(x_i)
        % if there is no end point, then it is a enclosed circle
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
    else
        [y_i, x_i]=find(end_points_i);
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
    end
    
    
    if(SaveFigures==1)
        h23 = figure(23);hold on;
        plot(line_i_x,line_i_y,'color',[100 100 255]/255);
        
          h26= figure(26);hold on;
    plot(line_i_x,line_i_y,'color',[100 100 255]/255);

    end
    
    
    
end

saveas(h23,[FilamentSegmentationChannelOutputDir,'/GEO/lightblue_all_', ...
    num2str(iFrame),'.fig']);


h24 = figure(24);hold off;
imagesc(imageInt_nohole);hold on; colormap(gray); axis image; hold on;



for i_area = changed_ind_good_I'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    % only for those good ones, check the curvature
    
    bw_i = zeros(size(bw_out));
    bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
    end_points_i = bwmorph(bw_i,'endpoints');
    [y_i, x_i]=find(end_points_i);
    
    if isempty(x_i)
        % if there is no end point, then it is a enclosed circle
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
    else
        [y_i, x_i]=find(end_points_i);
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
    end
    
    h24 = figure(24);hold on;
    plot(line_i_x,line_i_y,'c','LineWidth',1.1);
    
    
    
end


for i_area = changed_ind_good_K'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    % only for those good ones, check the curvature
    
    bw_i = zeros(size(bw_out));
    bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
    end_points_i = bwmorph(bw_i,'endpoints');
    [y_i, x_i]=find(end_points_i);
    
    if isempty(x_i)
        % if there is no end point, then it is a enclosed circle
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
    else
        [y_i, x_i]=find(end_points_i);
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
    end
    
    h24 = figure(24);hold on;
    plot(line_i_x,line_i_y,'color',[255 128 0   ] / 255,'LineWidth',1.1);
    
end
saveas(h24,[FilamentSegmentationChannelOutputDir,'/GEO/gm_all_', ...
    num2str(iFrame),'.fig']);

% get the mean intensity of the curves
for i_area = setdiff(setdiff(intersect(Bad_ind',ind_long'),changed_ind_good_K),changed_ind_good_I)
    [all_y_i, all_x_i] = find(labelMask == i_area);
    % only for those good ones, check the curvature
    
    bw_i = zeros(size(bw_out));
    bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
    end_points_i = bwmorph(bw_i,'endpoints');
    [y_i, x_i]=find(end_points_i);
    
    if isempty(x_i)
        % if there is no end point, then it is a enclosed circle
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
    else
        [y_i, x_i]=find(end_points_i);
        [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
    end
    
    h24 = figure(24);hold on;
    plot(line_i_x,line_i_y,'r');
    
        h26= figure(26);hold on;
    plot(line_i_x,line_i_y,'r');

end
axis off;
saveas(h24,[FilamentSegmentationChannelOutputDir,'/GEO/red_all_', ...
    num2str(iFrame),'.fig']);
saveas(h26,[FilamentSegmentationChannelOutputDir,'/GEO/light_bluered_all_', ...
    num2str(iFrame),'.fig']);
