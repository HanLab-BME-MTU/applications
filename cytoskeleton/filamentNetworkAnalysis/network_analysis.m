function [output_feature, VIF_ROI_model, VIF_ROI_orientation_model] = ...
    network_analysis(VIF_current_model,...
    VIF_current_seg, CellROI, radius,feature_flag)
% function for calculation the property of network, simply based the network structures.
% for input description see function load_MD_network_for_analysis.

% Liya Ding 01.2014.

% Initialize output feature struct
output_feature=[];

output_feature.straightness_per_filament_pool=[];
output_feature.curvature_per_pixel_pool=[];
output_feature.length_per_filament_pool=[];
output_feature.pixel_number_per_filament_pool=[];
output_feature.density_filament=[];
output_feature.scrabled_density_filament=[];
output_feature.orientation_pixel_pool_display=[];
output_feature.orientation_pixel_pool_display_center=[];


%% confine filament model in CellROI and also do digitalization on model
% size of image
img_size = size(VIF_current_seg);

% transfer filament network model to bw image
VIF_current_seg = filament_model_to_seg_bwim(VIF_current_model,img_size,[]);
              
[VIF_ROI_model, VIF_ROI_orientation_model] = ...
    ROIed_digital_filament_model(VIF_current_model,img_size, CellROI);

%% First check if needed to do analysis

% if not requested, nothing to do
if(sum(feature_flag(1:5))==0)
    return;
end

%% Start analysis

% define the initial pool array
orientation_pixel_pool = [];
pixel_number_per_filament_pool = [];
length_per_filament_pool = [];
straightness_per_filament_pool = [];
filament_mean_curvature = [];
curvature_per_pixel_pool =[];

%since these are calcuated fast and in same loop, calculate together

if(sum(feature_flag(1:3)>0) || feature_flag(6)>0 || feature_flag(17)>0 || feature_flag(18)>0 ) 
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
            
            
            if(feature_flag(17)>0 || feature_flag(18)>0)              
                
                line_smooth_H = fspecial('gaussian',5,1.5);
                
                line_i_x = (imfilter(x, line_smooth_H, 'replicate', 'same'));
                line_i_y = (imfilter(y, line_smooth_H, 'replicate', 'same'));
                
                Vertices = [line_i_x line_i_y];
                Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
                k=LineCurvature2D(Vertices,Lines);
                
                filament_mean_curvature(iF) = mean(abs(k));
                curvature_per_pixel_pool = [curvature_per_pixel_pool; k];
                
            end
            
            
            if(feature_flag(6)>0)
                orientation_pixel_pool = [orientation_pixel_pool; VIF_ROI_orientation_model{iF}];
            end
        end
    end
end

% Let the user deside later
% pixel_number_per_filament_pool = ...
%     pixel_number_per_filament_pool(pixel_number_per_filament_pool>=min_length);

% % length distribution, by the distance along the filament
% length_per_filament_pool = ...
%     length_per_filament_pool(length_per_filament_pool>=min_length);


% output according to user request
if(feature_flag(1)>0)
    
    output_feature.straightness_per_filament_pool=straightness_per_filament_pool;
end

if(feature_flag(2)>0)
    output_feature.length_per_filament_pool=length_per_filament_pool;
end

if(feature_flag(3)>0)
    output_feature.pixel_number_per_filament_pool=pixel_number_per_filament_pool;
end

%% density mesure
if(feature_flag(4)>0)    
    density_H = double((fspecial('disk',radius))>0);
    density_H = density_H./(sum(sum(density_H)));
    
    density_filament = imfilter(VIF_current_seg, density_H, 'replicate','same');
    T_sigma=40;
    O_sigma=pi/4;
    output_feature.density_filament=density_filament;
    output_feature.current_seg = VIF_current_seg;
end

%% scrambled density
if(feature_flag(5)>0)
    scrable_output_feature = ...
        scrable_network_analysis(VIF_current_model,VIF_current_seg,CellROI, radius,T_sigma,O_sigma);
    output_feature.scrabled_density_filament=scrable_output_feature.density_filament;
    output_feature.scrabled_filament=scrable_output_feature.current_seg;
end


if(feature_flag(6)>0)
    
    % Organize the orientation between -pi/2 ~ pi/2, with the corretion of pi/2
    
    orientation_pixel_pool_display = orientation_pixel_pool;
    orientation_pixel_pool_display = orientation_pixel_pool_display +pi/2;
    orientation_pixel_pool_display = mod(orientation_pixel_pool_display, pi);
    orientation_pixel_pool_display(orientation_pixel_pool_display>=pi/2) = orientation_pixel_pool_display(orientation_pixel_pool_display>=pi/2)-pi;
end


if(feature_flag(7)>0)
    
    % center the histogram at the mode
    
    [h,b] = hist(orientation_pixel_pool_display,-pi/2:pi/18:pi/2);
    h=h(1:end);
    h = h/length(orientation_pixel_pool_display);
    
    ind_max_h = find(h==max(h));
    ind_max_h = ind_max_h(1);
    mode_bin = b(ind_max_h)-pi/36;
    mode_bin = mod(mode_bin,pi);
    
    orientation_pixel_pool_display_center = orientation_pixel_pool_display-mode_bin;
    orientation_pixel_pool_display_center = mod(orientation_pixel_pool_display_center,pi);
    orientation_pixel_pool_display_center(orientation_pixel_pool_display_center>=pi/2) = orientation_pixel_pool_display_center(orientation_pixel_pool_display_center>=pi/2)-pi;
    
    % center the orientation at 0
%     orientation_pixel_pool_display_center = orientation_pixel_pool_display_center + mode_bin;
    orientation_pixel_pool_display_center = orientation_pixel_pool_display_center + 0;
else
    orientation_pixel_pool_display_center=[];
end
    
if(feature_flag(6)>0)
output_feature.orientation_pixel_pool_display=orientation_pixel_pool_display;
end

if(feature_flag(7)>0)
output_feature.orientation_pixel_pool_display_center=orientation_pixel_pool_display_center;
end

if(feature_flag(17)>0)
    output_feature.filament_mean_curvature =  filament_mean_curvature;   
end

if(feature_flag(18)>0)
    output_feature.curvature_per_pixel_pool = curvature_per_pixel_pool;
end


