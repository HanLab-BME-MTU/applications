function output_feature = perfilament_analysis(VIF_ROI_model, VIF_ROI_orientation_model,...
    VIF_current_seg, scale_map,current_img, MAX_st_res, feature_flag)
% function for calculation the property of network with other features like
% intensity, steerable filtering responce and etc.

% here the input need to be ROIed first, which is done in the
% network_analysis function.

% for input description see function load_MD_network_for_analysis.

% Liya Ding 2013

% initialization of output feature struct

output_feature=[];

output_feature.intensity_per_filament_pool = [];
output_feature.mean_intensity_per_filament_pool = [];
output_feature.intensity_per_fat_filament_pool = [];
output_feature.mean_intensity_per_fat_filament_pool = [];
output_feature.scale_per_filament_pool = [];

output_feature.st_per_filament_pool = [];
output_feature.mean_st_per_filament_pool = [];
output_feature.st_per_fat_filament_pool = [];
output_feature.mean_st_per_fat_filament_pool = [];

% if not requested, nothing to do
if(sum(feature_flag(8:end))==0)
    return;
end

%% start analysis
img_size = size(VIF_current_seg);

% Initialize features as empty ( since some are not going to be needed, so just empty)
intensity_per_filament_pool = [];
mean_intensity_per_filament_pool = [];
intensity_per_fat_filament_pool = [];
mean_intensity_per_fat_filament_pool = [];
scale_per_filament_pool = [];
st_per_filament_pool = [];
mean_st_per_filament_pool = [];
st_per_fat_filament_pool = [];
mean_st_per_fat_filament_pool = [];        

% for each filament, calculate the length, orientation, and straightness
% and put them into the pool
for iF = 1 : length(VIF_ROI_model)
    try
        x = VIF_ROI_model{iF}(:,1);
        y = VIF_ROI_model{iF}(:,2);
        
        filament_bw = zeros(img_size);
        filament_bw(sub2ind(img_size,round(y),round(x)))=1;
        fat_filament_bw = imdilate(filament_bw,ones(3,3));        
             
        if(sum(feature_flag(8:11))>0)
            % get the intenisty values
            filament_intensity = current_img(sub2ind(img_size,round(y),round(x)));
            fat_filament_intensity = current_img(fat_filament_bw>0);
            % gather intensity values
            intensity_per_filament_pool(iF) = sum(filament_intensity);
            mean_intensity_per_filament_pool(iF) = mean(filament_intensity);
            intensity_per_fat_filament_pool(iF) = sum(fat_filament_intensity);
            mean_intensity_per_fat_filament_pool(iF) = mean(fat_filament_intensity);
        end
        
        if(feature_flag(12)>0)
            scale_per_filament = scale_map(fat_filament_bw>0);
            scale_per_filament_pool(iF) = mean(scale_per_filament);
        end
        
        if(sum(feature_flag(13:16))>0)
            % get the steerable filtering responses
            filament_st = MAX_st_res(filament_bw>0);
            fat_filament_st = MAX_st_res(fat_filament_bw>0);
            
            % get the steerable filtering numbers
            st_per_filament_pool(iF) = sum(filament_st);
            mean_st_per_filament_pool(iF) = mean(filament_st);
            st_per_fat_filament_pool(iF) = sum(fat_filament_st);
            mean_st_per_fat_filament_pool(iF) = mean(fat_filament_st);
        end        
    end
end



% put the output into the feature structure, according to user request

if(feature_flag(8)>0)
output_feature.intensity_per_filament_pool = intensity_per_filament_pool;
end

if(feature_flag(9)>0)
output_feature.mean_intensity_per_filament_pool = mean_intensity_per_filament_pool;
end

if(feature_flag(10)>0)
output_feature.intensity_per_fat_filament_pool = intensity_per_fat_filament_pool;
end

if(feature_flag(11)>0)
output_feature.mean_intensity_per_fat_filament_pool = mean_intensity_per_fat_filament_pool;
end

if(feature_flag(12)>0)
output_feature.scale_per_filament_pool = scale_per_filament_pool;
end

if(feature_flag(13)>0)
output_feature.st_per_filament_pool = st_per_filament_pool;
end

if(feature_flag(14)>0)
output_feature.mean_st_per_filament_pool = mean_st_per_filament_pool;
end

if(feature_flag(15)>0)
output_feature.st_per_fat_filament_pool = st_per_fat_filament_pool;
end

if(feature_flag(16)>0)
output_feature.mean_st_per_fat_filament_pool = mean_st_per_fat_filament_pool;
end



