function output_feature = perfilament_analysis(VIF_current_model,...
    VIF_orientation, VIF_current_seg,outdir,iChannel,iFrame,...
    ROI, min_length,im_name,radius,scale_map,current_img, MAX_st_res)
% function for calculation the property of network with other features like
% intensity,steerable filtering responce and etc.

% Liya Ding 2013


% don't save every figure generated, unless debugging
save_everything_flag = 1;

img_size = size(VIF_current_seg);
if(iscell(im_name))
    im_name = im_name{1};
end

% current_img = imread(im_name);


im_name(im_name=='_')='-';
im_name_title{1} = im_name;

% transfer filament network model to bw image
VIF_current_seg = filament_model_to_seg_bwim(VIF_current_model,img_size,[]);

% transfer model to digital
[VIF_digital_model,VIF_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
    = filament_model_to_digital_with_orientation(VIF_current_model);


% if the ROI is full area, just copy
if(mean2(double(ROI))==1)
    VIF_ROI_model= VIF_digital_model;
    VIF_ROI_orientation_model =VIF_orientation_model;
else
    % if there is defined ROI, check filament by filament.
    count=1;
    VIF_ROI_model = cell(1,1);
    VIF_ROI_orientation_model =  cell(1,1);
    
    for iF = 1 : length(VIF_digital_model)
        try
            x = VIF_digital_model{iF}(:,1);
            y = VIF_digital_model{iF}(:,2);
            ind_xy = sub2ind(img_size,y,x);
            roi_flag_array = ROI(ind_xy);
            [inROI_label, inROI_N] = bwlabel(roi_flag_array);
            for iR = 1 : inROI_N
                if (sum(double(inROI_label==iR))>2)
                    VIF_ROI_model{count} = [x(find(inROI_label==iR)) y(find(inROI_label==iR))];
                    VIF_ROI_orientation_model{count} = VIF_orientation_model{iF}(find(inROI_label==iR));
                    count = count +1 ;
                end
            end
        end
    end
end

%  [Vif_ROIed_model,Vif_ROIed_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
%             = filament_model_to_digital_with_orientation(VIF_ROI_model);
%

% define the initial pool array
orientation_pixel_pool = [];
pixel_number_per_filament_pool = [];
length_per_filament_pool = [];
straightness_per_filament_pool = [];

% for each filament, calculate the length, orientation, and straightness
% and put them into the pool.
for iF = 1 : length(VIF_ROI_model)
    try
        x = VIF_ROI_model{iF}(:,1);
        y = VIF_ROI_model{iF}(:,2);
        
        filament_bw = zeros(img_size);
        filament_bw(sub2ind(img_size,round(y),round(x)))=1;
        fat_filament_bw = imdilate(filament_bw,ones(3,3));
        
        % get the intenisty values
        filament_intensity = current_img(sub2ind(img_size,round(y),round(x)));
        fat_filament_intensity = current_img(fat_filament_bw>0);
        
        % get the steerable filtering responses
        filament_st = MAX_st_res(filament_bw>0);
        fat_filament_st = MAX_st_res(fat_filament_bw>0);
        
        scale_per_filament = scale_map(fat_filament_bw>0);
        orientation_pixel_pool = [orientation_pixel_pool; VIF_ROI_orientation_model{iF}];
        
        % gather intensity values
        intensity_per_filament_pool(iF) = sum(filament_intensity);
        mean_intensity_per_filament_pool(iF) = mean(filament_intensity);
        intensity_per_fat_filament_pool(iF) = sum(fat_filament_intensity);
        mean_intensity_per_fat_filament_pool(iF) = mean(fat_filament_intensity);
        
        % get the steerable filtering numbers
        st_per_filament_pool(iF) = sum(filament_st);
        mean_st_per_filament_pool(iF) = mean(filament_st);
        st_per_fat_filament_pool(iF) = sum(fat_filament_st);
        mean_st_per_fat_filament_pool(iF) = mean(fat_filament_st);
        
        scale_per_filament_pool(iF) = mean(scale_per_filament);
        
    end
end

pixel_number_per_filament_pool = ...
    pixel_number_per_filament_pool(pixel_number_per_filament_pool>=min_length);


% length distribution, by the distance along the filament
length_per_filament_pool = ...
    length_per_filament_pool(length_per_filament_pool>=min_length);

% Organize the orientation between -pi/2 ~ pi/2, with the corretion of pi/2

orientation_pixel_pool_display = orientation_pixel_pool;
orientation_pixel_pool_display = orientation_pixel_pool_display +pi/2;
orientation_pixel_pool_display = mod(orientation_pixel_pool_display, pi);
orientation_pixel_pool_display(orientation_pixel_pool_display>=pi/2) = orientation_pixel_pool_display(orientation_pixel_pool_display>=pi/2)-pi;

figure_flag=1;
if(figure_flag>0)
    h3 =  figure(3);
    rose(orientation_pixel_pool_display);
    title([im_name_title,' Orientation of Filaments']);
    
    saveas(h3, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h3, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
    saveas(h3, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
    
    
    % the orientation in histogram
    h6 =  figure(6);
    
    [h,b] = hist(orientation_pixel_pool_display,-pi/2:pi/18:pi/2);
    h=h(1:end);
    h = h/length(orientation_pixel_pool_display);
    bin = (-pi/2:pi/18:pi/2) +pi/36;
    bar(bin,h);
    axis([-pi/2-pi/36 pi/2+pi/36 0 0.3]);
    
    title([im_name_title,' Orientation of Filaments']);
    
    saveas(h6, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h6, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
    saveas(h6, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
    
    
% center the histogram at the mode    
    ind_max_h = find(h==max(h));
ind_max_h = ind_max_h(1);
mode_bin = b(ind_max_h)-pi/36;
    mode_bin = mod(mode_bin,pi);
    
orientation_pixel_pool_display_center = orientation_pixel_pool_display-mode_bin;
orientation_pixel_pool_display_center = mod(orientation_pixel_pool_display_center,pi);
orientation_pixel_pool_display_center(orientation_pixel_pool_display_center>=pi/2) = orientation_pixel_pool_display_center(orientation_pixel_pool_display_center>=pi/2)-pi;

orientation_pixel_pool_display_center = orientation_pixel_pool_display_center + mode_bin;


 h13 =  figure(13);
    rose(orientation_pixel_pool_display_center);
    title([im_name_title,' Orientation of Filaments']);
    
    saveas(h13, [outdir,filesep,'network_orientationrose_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h13, [outdir,filesep,'network_orientationrose_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
   
    
    % the orientation in histogram
    h16 =  figure(16);
    
    [h,b] = hist(orientation_pixel_pool_display_center,-pi/2+mode_bin:pi/18:pi/2+mode_bin);
    h=h(1:end);
    h = h/length(orientation_pixel_pool_display_center);
    bin = (-pi/2:pi/18:pi/2)+mode_bin+pi/36;
    bar(bin,h);
    axis([-pi/2-pi/36+mode_bin pi/2+pi/36+mode_bin 0 0.3]);
    
    title([im_name_title,' Orientation of Filaments']);
    
    saveas(h16, [outdir,filesep,'network_orientationhist_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h16, [outdir,filesep,'network_orientationhist_centered_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
   


end

% put the output into the feature structure
output_feature=[];

output_feature.orientation_pixel_pool_display=orientation_pixel_pool_display;
output_feature.orientation_pixel_pool_display_center=orientation_pixel_pool_display_center;
output_feature.intensity_per_filament_pool = intensity_per_filament_pool;
output_feature.mean_intensity_per_filament_pool = mean_intensity_per_filament_pool;
output_feature.intensity_per_fat_filament_pool = intensity_per_fat_filament_pool;
output_feature.mean_intensity_per_fat_filament_pool = mean_intensity_per_fat_filament_pool;
output_feature.scale_per_filament_pool = scale_per_filament_pool;

output_feature.st_per_filament_pool = st_per_filament_pool;
output_feature.mean_st_per_filament_pool = mean_st_per_filament_pool;
output_feature.st_per_fat_filament_pool = st_per_fat_filament_pool;
output_feature.mean_st_per_fat_filament_pool = mean_st_per_fat_filament_pool;


