function output_feature = network_analysis(VIF_current_model,...
    VIF_orientation, VIF_current_seg,outdir,iChannel,iFrame,...
    ROI, min_length,im_name,radius)
% function for calculation the property of network, simply based the network structures.

% Liya Ding 01.2014.

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
        
        filament_detailed_length = sum(sqrt((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2));
        filament_start_end_distance = sqrt((x(1)-x(end)).^2 + (y(1)-y(end)).^2);
        
        pixel_number_per_filament_pool(iF) = length(x);
        length_per_filament_pool(iF) = filament_detailed_length;
        straightness_per_filament_pool(iF) = filament_start_end_distance/filament_detailed_length;
        
    end
end

pixel_number_per_filament_pool = ...
    pixel_number_per_filament_pool(pixel_number_per_filament_pool>=min_length);


% Plot these pools out
figure_flag=1;

if(figure_flag>0)
    % length distribution, by how many pixels
    h1 =  figure(1);
    
    [h,bin] = hist(pixel_number_per_filament_pool,0:20:1000);
    h = h/length(pixel_number_per_filament_pool);
    bar(bin,h);
    axis([-10 310 0 0.3]);
    
    title([im_name_title,' Pixels Number Distribution']);
    saveas(h1, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h1, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
    saveas(h1, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
end

% length distribution, by the distance along the filament
length_per_filament_pool = ...
    length_per_filament_pool(length_per_filament_pool>=min_length);
[h,bin] = hist(length_per_filament_pool,0:20:1000);
h = h/length(length_per_filament_pool);

if(figure_flag>0)
    
    h2 =  figure(2);
    
    
    bar(bin,h);
    axis([-10 310 0 0.3]);
    
    title([im_name_title,' Length Distribution']);
    saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
    saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
    
    
    % the straightness distribution
    h4 =  figure(4);
    [h,bin] = hist(straightness_per_filament_pool,0:0.02:1);
    h = h/length(straightness_per_filament_pool);
    bar(bin,h);
    axis([0.69 0.97 0 0.2]);
    
    title([im_name_title,' Straightness']);
    
    saveas(h4, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h4, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
    saveas(h4, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
    
    h5 =  figure(5);
    h = boxplot(straightness_per_filament_pool);
    set(h(7,:),'Visible','off');
    axis([0 2 0.6 1.01]);
    
    title([im_name_title,' Straightness']);
    
    saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
    saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
    saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);
end

density_H = double((fspecial('disk',radius))>0);
density_H = density_H./(sum(sum(density_H)));

density_filament = imfilter(VIF_current_seg,density_H, 'replicate','same');
T_sigma=40;
O_sigma=pi/4;
scrable_output_feature = ...
    scrable_network_analysis(VIF_current_model,VIF_current_seg,ROI, radius,T_sigma,O_sigma);


output_feature=[];

output_feature.straightness_per_filament_pool=straightness_per_filament_pool;
output_feature.length_per_filament_pool=length_per_filament_pool;
output_feature.pixel_number_per_filament_pool=pixel_number_per_filament_pool;
output_feature.density_filament=density_filament;
output_feature.scrabled_density_filament=scrable_output_feature.density_filament;




