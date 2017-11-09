% % function to do branch analysis for locally copied folders fo branch
% % analysis results
% % Liya Ding, March, 2015
% 
 nList = 4;
% BA_Group = cell(1, nList);

T_branchsize=200;
T_branchduration=2;
Group_ROOT_DIR=[];
branch_percent_threshold = [ 5 10 20]; % anything larger than 1% of cell size
short_time_frames = [1 3 5]; % check things on each period of time  1 frame
speed_percent_threshold = [5 10]; % the speed to be considered as valid, not reorienting, 1% of cell demension( sqrt of cell pixels) is very low
branch_orientation_p_threshold = [ 0.999]; %  the threshold to consider branch orientation as different 0.99 means they can be really similar to each other
T_speed=4;
ROOT_DIR = 'F:\Work\fromNancy\branch_analysis\branchsize_framebyframe';
% 
% close all;
% 
% for iML = 1 : nList    
%  
%     MD_name_array = folder_name_list{iML};
%     folder_name = MD_name_array{1};  
%     
%     ind__s = find(folder_name=='_');
%     ML_name = folder_name(1:ind__s(2)-1);
%     
%     ML_ROOT_DIR = [ROOT_DIR, filesep, ML_name,'_plot'];
%     
%     % the number of movies
%     cellNumber =  numel(folder_name_list{iML});
%     
%     BA_output_ML_cell= cell(1,1,1);
%     
%     for iCell = 1: cellNumber
%         display(['Checking: iCell:', num2str(iCell)]);
%         
%         % the folder name if there is marking
%         outputPath = [ROOT_DIR,filesep,MD_name_array{iCell}];
%         
%         % check if this folder exist
%         if(exist(outputPath,'dir'))
%             % if it exist, try to do the branch analysis
%             try
%                 if(exist([outputPath,filesep,'branch_analysis_results_balloon.mat'],'file'))
%                     load([outputPath,filesep,'branch_analysis_results_balloon.mat'],'BA_output');
%                 else
%                     continue;
%                 end
%                 
%                 folder_name = MD_name_array{iCell};
%                 ind__s = find(folder_name=='_');
%                 Identifier = folder_name(1:ind__s(3)-1);
%                 BA_output.Identifier = Identifier;
%                 % save
%                 BA_output_ML_cell{1, 1}{1, iCell} = BA_output;
%             end
%         end
%         BA_Group{1, iML} = BA_output_ML_cell;
%     end
% end


%%

% initialize the pools

% Group_Pool_Travel_Length = [];
% Group_Pool_Travel_Distance = [];
% Group_Pool_Travel_Speed = [];
% Group_Pool_Travel_Speed_without_pausing = [];
% Group_Pool_Cell_Marked_Frame_Number = [];
% Group_Pool_branch_number_tracked = [];
% Group_Pool_branch_number_max = [];
% Group_Pool_branch_number_mean = [];
% Group_Pool_branch_duration_array = [];
% Group_Pool_branch_duration_mean= [];
% Group_Pool_branch_vif_mean_intensity= [];
% Group_Pool_protrusion_vif_mean_intensity= [];
% Group_Pool_retraction_vif_mean_intensity= [];
% Group_Pool_whole_cell_vif_mean_intensity= [];
% Group_Pool_whole_cell_size_mean= [];
% Group_Pool_branch_number_mean_pat=[];
% Group_Pool_whole_cell_vim_totalamount_mean=[];
% Group_Pool_whole_cell_vif_mean_intensity_pat=[];
% Group_Pool_branch_size_mean=[];
% Group_Pool_thresholded_branch_number_mean=[];
% Group_Pool_branch_cellmovement_std=[];
% Group_Pool_fila_branch_orientation_pool_std=[];
% Group_Pool_fila_trajectory_orientation_pool_std=[];
% Group_Pool_fila_branch_orientation_pool=[];
% Group_Pool_fila_trajectory_orientation_pool=[];
% Group_Pool_fila_branch_orientation_pool_slow=[];
% Group_Pool_fila_trajectory_orientation_pool_slow=[];
% Group_Pool_fila_branch_orientation_pool_fast=[];
% Group_Pool_fila_trajectory_orientation_pool_fast=[];
% Group_Pool_pool_vif_intensity=[];
% Group_Pool_CompletedFrame_last=[];
% Group_Pool_whole_cell_vim_seg_total = [];
% Group_Pool_whole_cell_vim_seg_mean  = [];
% Group_Pool_whole_cell_vim_nms_total = [];
% Group_Pool_whole_cell_vim_nms_mean  = [];
% Group_Pool_branch_seg_total = [];
% Group_Pool_branch_seg_mean  = [];
% Group_Pool_branch_nms_total = [];
% Group_Pool_branch_nms_mean  = [];
% Group_Pool_branch_number_weighted = cell(1,numel(branch_percent_threshold));
% Group_Pool_Travel_Speed_short_time = cell(1, numel(short_time_frames));
% Group_Pool_reorienting_status_short_time = cell(numel(short_time_frames),numel(speed_percent_threshold),numel(branch_orientation_p_threshold));
% Group_Pool_vim_fila_shorttime = cell(1, numel(short_time_frames));
% Group_Pool_vim_nms_shorttime = cell(1, numel(short_time_frames));
% Group_Pool_vim_int_shorttime = cell(1, numel(short_time_frames));
% Group_Pool_persistency_shorttime  = cell(1, numel(short_time_frames));

Group_Pool_curvature_map_pool = [];
Group_Pool_filament_alignment_map_pool = [];
Group_Pool_filament_density_pool = [];
Group_Pool_prot_or_retr_fila_pool = [];
Group_Pool_curvature_map_full_pool = [];
Group_Pool_filament_alignment_map_full_pool = [];
Group_Pool_filament_density_full_pool = [];
Group_Pool_prot_or_retr_fila_full_pool = [];
Group_Pool_scale_map_pool = [];

Group_Pool_last_curvature_map_pool = [];
Group_Pool_last_filament_alignment_map_pool = [];
Group_Pool_last_filament_density_pool = [];
Group_Pool_last_prot_or_retr_fila_pool = [];
Group_Pool_last_curvature_map_full_pool = [];
Group_Pool_last_filament_alignment_map_full_pool = [];
Group_Pool_last_filament_density_full_pool = [];
Group_Pool_last_prot_or_retr_fila_full_pool = [];
Group_Pool_last_scale_map_pool = [];

Group_Pool_int_map_pool = [];
Group_Pool_st_map_pool = [];
Group_Pool_int_map_full_pool = [];
Group_Pool_st_map_full_pool = [];

Group_Pool_last_int_map_pool = [];
Group_Pool_last_st_map_pool = [];
Group_Pool_last_int_map_full_pool = [];
Group_Pool_last_st_map_full_pool = [];


Identifier_cell = [];

iAllCell = 0;

% if no input, the DIR of last movieList
Group_ROOT_DIR = [ROOT_DIR,filesep,'All_4_days'];

if(exist(Group_ROOT_DIR,'dir')==0)
    mkdir(Group_ROOT_DIR);
end

%% For each ML
for iML = 1 : 4
    BA_output_ML_cell = BA_Group{1, iML};
    
%     ML_Pool_Travel_Length = [];
%     ML_Pool_Travel_Distance = [];
%     ML_Pool_Travel_Speed = [];
%     ML_Pool_Cell_Marked_Frame_Number = [];
%     ML_Pool_branch_number_tracked = [];
%     ML_Pool_branch_number_max = [];
%     ML_Pool_branch_number_mean = [];
%     ML_Pool_branch_duration_array = [];
%     ML_Pool_branch_duration_mean= [];
%     ML_Pool_branch_vif_mean_intensity= [];
%     ML_Pool_protrusion_vif_mean_intensity= [];
%     ML_Pool_retraction_vif_mean_intensity= [];
%     ML_Pool_whole_cell_vif_mean_intensity= [];
%     ML_Pool_whole_cell_size_mean= [];
%     ML_Pool_branch_number_mean_pat=[];
%     ML_Pool_whole_cell_vim_totalamount_mean=[];
%     ML_Pool_whole_cell_vif_mean_intensity_pat=[];
%     ML_Pool_branch_size_mean=[];
%     ML_Pool_thresholded_branch_number_mean=[];
%     ML_Pool_branch_cellmovement_std=[];
%     ML_Pool_fila_branch_orientation_pool_std=[];
%     ML_Pool_fila_trajectory_orientation_pool_std=[];
%     ML_Pool_fila_branch_orientation_pool=[];
%     ML_Pool_fila_trajectory_orientation_pool=[];
%     ML_Pool_fila_branch_orientation_pool_slow=[];
%     ML_Pool_fila_trajectory_orientation_pool_slow=[];
%     ML_Pool_fila_branch_orientation_pool_fast=[];
%     ML_Pool_fila_trajectory_orientation_pool_fast=[];
%     ML_Pool_CompletedFrame_last = [];
%     ML_Pool_pool_vif_intensity=[];
%     ML_Pool_whole_cell_vim_seg_total = [];
%     ML_Pool_whole_cell_vim_seg_mean = [];
%     ML_Pool_whole_cell_vim_nms_total = [];
%     ML_Pool_whole_cell_vim_nms_mean = [];
%     ML_Pool_branch_seg_total = [];
%     ML_Pool_branch_seg_mean = [];
%     ML_Pool_branch_nms_total = [];
%     ML_Pool_branch_nms_mean = [];
%     ML_Pool_Travel_Speed_without_pausing = [];
%     
%     ML_Pool_branch_number_weighted = cell(1,numel(branch_percent_threshold));
%     ML_Pool_Travel_Speed_short_time = cell(1, numel(short_time_frames));
%     ML_Pool_reorienting_status_short_time = cell(numel(short_time_frames),numel(speed_percent_threshold),numel(branch_orientation_p_threshold));
%     ML_Pool_vim_fila_shorttime = cell(1, numel(short_time_frames));
%     ML_Pool_vim_nms_shorttime = cell(1, numel(short_time_frames));
%     ML_Pool_vim_int_shorttime = cell(1, numel(short_time_frames));
%     ML_Pool_persistency_shorttime  = cell(1, numel(short_time_frames));
%     
    
     ML_Pool_curvature_map_pool = [];
    ML_Pool_filament_alignment_map_pool = [];
    ML_Pool_filament_density_pool = [];
    ML_Pool_prot_or_retr_fila_pool = [];
    ML_Pool_curvature_map_full_pool = [];
    ML_Pool_filament_alignment_map_full_pool = [];
    ML_Pool_filament_density_full_pool = [];
    ML_Pool_prot_or_retr_fila_full_pool = [];
    ML_Pool_scale_map_pool = [];
    
    ML_Pool_last_curvature_map_pool = [];
    ML_Pool_last_filament_alignment_map_pool = [];
    ML_Pool_last_filament_density_pool = [];
    ML_Pool_last_prot_or_retr_fila_pool = [];
    ML_Pool_last_curvature_map_full_pool = [];
    ML_Pool_last_filament_alignment_map_full_pool = [];
    ML_Pool_last_filament_density_full_pool = [];
    ML_Pool_last_prot_or_retr_fila_full_pool = [];
    ML_Pool_last_scale_map_pool = [];
    
    ML_Pool_int_map_pool = [];
    ML_Pool_st_map_pool = [];
    ML_Pool_int_map_full_pool = [];
    ML_Pool_st_map_full_pool = [];
    
    ML_Pool_last_int_map_pool = [];
    ML_Pool_last_st_map_pool = [];
    ML_Pool_last_int_map_full_pool = [];
    ML_Pool_last_st_map_full_pool = [];
    
    
    MD_name_array = folder_name_list{iML};
    folder_name = MD_name_array{1};
    ind__s = find(folder_name=='_');
    ML_name = folder_name(1:ind__s(2)-1);
    
    ML_ROOT_DIR = [ROOT_DIR, filesep, ML_name,'_plot'];
    
    % the number of movies
    cellNumber =  numel(MD_name_array);
    
    for iCell = 1 :  cellNumber
        try
            BA_output = BA_output_ML_cell{1, 1}{1, iCell};
        catch
            continue;
        end
        
        if(~isempty(BA_output))
%             if isfield(BA_output,'speed_marked_frames') && isfield(BA_output,'cell_vif_seg_total_array') ...
%                     && numel(BA_output.cell_vif_nms_total_array' )== numel(BA_output.cell_size_array)
%                 S = BA_output.speed_marked_frames;
%                 ML_Pool_Travel_Speed_without_pausing = [ML_Pool_Travel_Speed_without_pausing mean(S(S>T_speed))];
%                 smoothed_speed_marked_frames = BA_output.smoothed_speed_marked_frames;
%                 center_x = BA_output.center_x;
%                 center_y = BA_output.center_y;
%                 center_x = center_x(:);
%                 center_y = center_y(:);
%                 % assume the last frame is moving in the same speed
%                 % as the end-1 to end
%                 smoothed_speed_marked_frames(numel(center_x)) = smoothed_speed_marked_frames(numel(center_x)-1);
%                 smoothed_speed_marked_frames = smoothed_speed_marked_frames(:);
%                 
%                 ML_Pool_pool_vif_intensity= [ML_Pool_pool_vif_intensity; BA_output.pool_all_vif_intensity];
%                    
% %                 try
%                 %%  speed with shorttime definition
%                 for iShorttime = 1 : numel(short_time_frames)
%                     ShortTimeFrames = short_time_frames(iShorttime);
%                     
%                     half_shorttime_frame = floor(ShortTimeFrames/2);
%                     smooth_center_x = BA_output.smooth_center_x;
%                     smooth_center_y = BA_output.smooth_center_y;
%                     
%                     persistency_shorttime=[];
%                     if(ShortTimeFrames>1)
%                         for iFrame = 1 : numel(smooth_center_x)
%                             
%                             this_x = smooth_center_x(max(1,iFrame-half_shorttime_frame):min(iFrame+half_shorttime_frame,numel(smooth_center_x)));
%                             this_y = smooth_center_y(max(1,iFrame-half_shorttime_frame):min(iFrame+half_shorttime_frame,numel(smooth_center_x)));
%                             if(numel(this_x)>1)
%                             direct_distance = sqrt((this_x(end)-this_x(1)).^2+(this_y(end)-this_y(1)).^2);
%                             path_distance = sum(sqrt((this_x(2:end)-this_x(1:end-1)).^2+(this_y(2:end)-this_y(1:end-1)).^2));
%                             
%                             persistency_shorttime(iFrame,1)=direct_distance/path_distance;
%                             else
%                                 persistency_shorttime(iFrame,1)=nan;
%                             end
%                         end
%                     else
%                         persistency_shorttime = nan(numel(smooth_center_x),1);
%                     end
%                     
%                     half_shorttime_frame = 1;
%                     
%                     %                             half_shorttime_frame = floor(ShortTimeFrames/2);
%                    
%                     H = ones(ShortTimeFrames,1)/ShortTimeFrames;
%                     
%                     smoothed_speed_marked_frames_shorttime = imfilter(smoothed_speed_marked_frames,H,'same','replicate');
%                     
%                     
%                     vif_fila_mean = BA_output.cell_vif_seg_total_array' ./(BA_output.cell_size_array);
%                     vif_nms_mean = BA_output.cell_vif_nms_total_array' ./(BA_output.cell_size_array);
%                     vif_int_mean = BA_output.cell_vimtotal_pool ./(BA_output.cell_size_array);
%                     
%                     vif_fila_mean = vif_fila_mean(:);
%                     vif_nms_mean  = vif_nms_mean(:);
%                     vif_int_mean  = vif_int_mean(:);
%                     
%                     vif_fila_mean_shorttime = imfilter(vif_fila_mean,H,'same','replicate');
%                     vif_nms_mean_shorttime = imfilter(vif_nms_mean,H,'same','replicate');
%                     vif_int_mean_shorttime = imfilter(vif_int_mean,H,'same','replicate');
%                     
%                     ML_Pool_vim_fila_shorttime{iShorttime} = [ML_Pool_vim_fila_shorttime{iShorttime}; vif_fila_mean_shorttime(half_shorttime_frame:end-half_shorttime_frame+1,1)];
%                     ML_Pool_vim_nms_shorttime{iShorttime} = [ML_Pool_vim_nms_shorttime{iShorttime}; vif_nms_mean_shorttime(half_shorttime_frame:end-half_shorttime_frame+1,1)];
%                     ML_Pool_vim_int_shorttime{iShorttime} = [ML_Pool_vim_int_shorttime{iShorttime}; vif_int_mean_shorttime(half_shorttime_frame:end-half_shorttime_frame+1,1)];
%                     ML_Pool_Travel_Speed_short_time{iShorttime} = [ML_Pool_Travel_Speed_short_time{iShorttime}; smoothed_speed_marked_frames_shorttime(half_shorttime_frame:end-half_shorttime_frame+1,1)] ;
%                     ML_Pool_persistency_shorttime{iShorttime} = [ML_Pool_persistency_shorttime{iShorttime}; persistency_shorttime(half_shorttime_frame:end-half_shorttime_frame+1,1)];
%                     
%                 end
%                 
%                 
%             ML_Pool_CompletedFrame_last=[ML_Pool_CompletedFrame_last BA_output.CompletedFrame(end)];
%             ML_Pool_branch_size_mean = [ML_Pool_branch_size_mean BA_output.branch_mean_size];
%             ML_Pool_whole_cell_size_mean= [ML_Pool_whole_cell_size_mean BA_output.whole_cell_size_mean];
%             ML_Pool_branch_number_mean_pat = [ML_Pool_branch_number_mean_pat repmat(BA_output.branch_number_mean, [1, length(BA_output.branch_vif_mean_intensity)])];
%             Identifier_cell{1,numel(Identifier_cell)+1} = BA_output.Identifier;
%             
%             Thresholded_branch_number_mean = sum(BA_output.branch_duration_array(find(BA_output.branch_mean_size>T_branchsize & BA_output.branch_duration_array>T_branchduration)))...
%                 ./ BA_output.cell_marked_frame_number;
%             ML_Pool_thresholded_branch_number_mean = [ML_Pool_thresholded_branch_number_mean Thresholded_branch_number_mean];
%             
%             
%                 
%                 
%                 %% branch number with weight
%                 for iB = 1 : numel(branch_percent_threshold)
%                     branch_size_percent = BA_output.branch_size_matrix./(repmat(BA_output.cell_size_array, 1, size(BA_output.branch_size_matrix,2)))*100;
%                     branch_size_percent(branch_size_percent> branch_percent_threshold(iB)) = branch_percent_threshold(iB);
%                     branch_size_percent = branch_size_percent/branch_size_percent(iB);
%                     branch_number_afterweight = sum(branch_size_percent,2);
%                     ML_Pool_branch_number_weighted{iB} = [ ML_Pool_branch_number_weighted{iB}; branch_number_afterweight];
%                 end
%                 %%
%                 
%                 for iSpeed = 1 : numel(speed_percent_threshold)
%                     for iPvalue = 1 : numel(branch_orientation_p_threshold)
%                         for iShorttime = 1 : numel(short_time_frames)
%                             ShortTimeFrames = short_time_frames(iShorttime);
%                                                         
%                             H = ones(ShortTimeFrames,1)/ShortTimeFrames;
%                             smoothed_speed_marked_frames_shorttime = imfilter(smoothed_speed_marked_frames,H,'same','replicate');
%                             
%                              this_speed = smoothed_speed_marked_frames_shorttime;
%                              this_speed_percent = this_speed./(sqrt(BA_output.cell_size_array))*100;
% %                              pvalue_array(1)=1;
% %                              for iFrame = 2 : numel(BA_output.branch_orienation_perframe_allbranch)
% %                              [Hipo, pValue, KSstatistic] = kstest2(BA_output.branch_orienation_perframe_allbranch{iFrame},BA_output.branch_orienation_perframe_allbranch{iFrame-1});
% %                              pvalue_array(iFrame) = pValue;
% %                              end
% %                              pvalue_array = pvalue_array(:);
%                              
%                              this_status = (this_speed_percent<speed_percent_threshold(iSpeed));
%                              
%                              ML_Pool_reorienting_status_short_time{iShorttime,iSpeed,iPvalue} = ...
%                                  [ML_Pool_reorienting_status_short_time{iShorttime,iSpeed,iPvalue};...
%                                  this_status];
%    
%                         end
%                     end
%                 end  
% %             end
%             end
%             
%             
            
                                
                    ML_Pool_curvature_map_pool= [ML_Pool_curvature_map_pool; BA_output.curvature_map_pool];
                    ML_Pool_filament_alignment_map_pool= [ML_Pool_filament_alignment_map_pool; BA_output.filament_alignment_map_pool];
                    ML_Pool_filament_density_pool= [ML_Pool_filament_density_pool; BA_output.filament_density_pool];
                    ML_Pool_prot_or_retr_fila_pool= [ML_Pool_prot_or_retr_fila_pool; BA_output.prot_or_retr_fila_pool];
                    ML_Pool_curvature_map_full_pool= [ML_Pool_curvature_map_full_pool; BA_output.curvature_map_full_pool];
                    ML_Pool_filament_alignment_map_full_pool= [ML_Pool_filament_alignment_map_full_pool; BA_output.filament_alignment_map_full_pool];
                    
                    ML_Pool_prot_or_retr_fila_full_pool= [ML_Pool_prot_or_retr_fila_full_pool; BA_output.prot_or_retr_fila_full_pool];
                    ML_Pool_scale_map_pool= [ML_Pool_scale_map_pool; BA_output.scale_map_pool];
                    
                    ML_Pool_last_curvature_map_pool= [ML_Pool_last_curvature_map_pool; BA_output.last_curvature_map_pool];
                    ML_Pool_last_filament_alignment_map_pool= [ML_Pool_last_filament_alignment_map_pool; BA_output.last_filament_alignment_map_pool];
                    ML_Pool_last_filament_density_pool= [ML_Pool_last_filament_density_pool; BA_output.last_filament_density_pool];
                    ML_Pool_last_prot_or_retr_fila_pool= [ML_Pool_last_prot_or_retr_fila_pool; BA_output.last_prot_or_retr_fila_pool];
                    ML_Pool_last_curvature_map_full_pool= [ML_Pool_last_curvature_map_full_pool; BA_output.last_curvature_map_full_pool];
                    ML_Pool_last_filament_alignment_map_full_pool= [ML_Pool_last_filament_alignment_map_full_pool; BA_output.last_filament_alignment_map_full_pool];
                    ML_Pool_last_filament_density_full_pool= [ML_Pool_last_filament_density_full_pool; BA_output.last_filament_density_full_pool];
                    ML_Pool_last_scale_map_pool= [ML_Pool_last_scale_map_pool; BA_output.last_scale_map_pool];
                    
                    
                    
                    ML_Pool_int_map_pool= [ML_Pool_int_map_pool; BA_output.int_map_pool];
                    ML_Pool_st_map_pool= [ML_Pool_st_map_pool; BA_output.st_map_pool];
                    ML_Pool_int_map_full_pool= [ML_Pool_int_map_full_pool; BA_output.int_map_full_pool];
                    ML_Pool_st_map_full_pool= [ML_Pool_st_map_full_pool; BA_output.st_map_full_pool];
                    
                    ML_Pool_last_int_map_pool= [ML_Pool_last_int_map_pool; BA_output.last_int_map_pool];
                    ML_Pool_last_st_map_pool= [ML_Pool_last_st_map_pool; BA_output.last_st_map_pool];
                    ML_Pool_last_int_map_full_pool= [ML_Pool_last_int_map_full_pool; BA_output.last_int_map_full_pool];
                    ML_Pool_last_st_map_full_pool= [ML_Pool_last_st_map_full_pool; BA_output.last_st_map_full_pool];
                    
                    
        
            iAllCell = iAllCell + 1;
        end
    end



%% vim normalization
vim_max = prctile(double(ML_Pool_pool_vif_intensity),98);
% 
% display('50')
% prctile(double(ML_Pool_pool_vif_intensity),50)
% display('75')
% prctile(double(ML_Pool_pool_vif_intensity),75)
% display('95')
% prctile(double(ML_Pool_pool_vif_intensity),95)
% display('98')
% prctile(double(ML_Pool_pool_vif_intensity),98)
% display('99')
% prctile(double(ML_Pool_pool_vif_intensity),99)
% 
% 
% ML_Pool_branch_vif_mean_intensity = ML_Pool_branch_vif_mean_intensity/vim_max;
% ML_Pool_protrusion_vif_mean_intensity = ML_Pool_protrusion_vif_mean_intensity/vim_max;
% ML_Pool_retraction_vif_mean_intensity = ML_Pool_retraction_vif_mean_intensity/vim_max;
% ML_Pool_whole_cell_vif_mean_intensity = ML_Pool_whole_cell_vif_mean_intensity/vim_max;
% ML_Pool_whole_cell_vif_mean_intensity_pat = ML_Pool_whole_cell_vif_mean_intensity_pat/vim_max;
% ML_Pool_whole_cell_vim_totalamount_mean = ML_Pool_whole_cell_vim_totalamount_mean/vim_max;
% for iShorttime = 1 : numel(short_time_frames)
%     ML_Pool_vim_int_shorttime{iShorttime} = ML_Pool_vim_int_shorttime{iShorttime}/vim_max;
% end

save([ML_name,'ML_Pools.mat']);
    
    %% in the end put data from this ML to the pool of all MLs
    
%     Group_Pool_Travel_Length = [Group_Pool_Travel_Length ML_Pool_Travel_Length];
%     Group_Pool_Travel_Distance = [Group_Pool_Travel_Distance ML_Pool_Travel_Distance];
%     Group_Pool_Travel_Speed = [Group_Pool_Travel_Speed ML_Pool_Travel_Speed];
%     Group_Pool_Travel_Speed_without_pausing = [Group_Pool_Travel_Speed_without_pausing ML_Pool_Travel_Speed_without_pausing];
%     Group_Pool_Cell_Marked_Frame_Number = [Group_Pool_Cell_Marked_Frame_Number ML_Pool_Cell_Marked_Frame_Number];
%     Group_Pool_branch_number_tracked = [Group_Pool_branch_number_tracked ML_Pool_branch_number_tracked];
%     Group_Pool_branch_number_max = [Group_Pool_branch_number_max ML_Pool_branch_number_max];
%     Group_Pool_branch_number_mean = [Group_Pool_branch_number_mean ML_Pool_branch_number_mean];
%     Group_Pool_branch_duration_array = [Group_Pool_branch_duration_array ML_Pool_branch_duration_array];
%     Group_Pool_branch_duration_mean= [Group_Pool_branch_duration_mean ML_Pool_branch_duration_mean];
%     Group_Pool_branch_vif_mean_intensity= [Group_Pool_branch_vif_mean_intensity ML_Pool_branch_vif_mean_intensity];
%     Group_Pool_protrusion_vif_mean_intensity= [Group_Pool_protrusion_vif_mean_intensity ML_Pool_protrusion_vif_mean_intensity];
%     Group_Pool_retraction_vif_mean_intensity= [Group_Pool_retraction_vif_mean_intensity ML_Pool_retraction_vif_mean_intensity];
%     Group_Pool_whole_cell_vif_mean_intensity= [Group_Pool_whole_cell_vif_mean_intensity ML_Pool_whole_cell_vif_mean_intensity];
%     Group_Pool_whole_cell_size_mean= [Group_Pool_whole_cell_size_mean ML_Pool_whole_cell_size_mean];
%     Group_Pool_branch_number_mean_pat=[Group_Pool_branch_number_mean_pat ML_Pool_branch_number_mean_pat];
%     Group_Pool_whole_cell_vim_totalamount_mean=[Group_Pool_whole_cell_vim_totalamount_mean ML_Pool_whole_cell_vim_totalamount_mean];
%     Group_Pool_whole_cell_vif_mean_intensity_pat=[Group_Pool_whole_cell_vif_mean_intensity_pat ML_Pool_whole_cell_vif_mean_intensity_pat];
%     Group_Pool_branch_size_mean=[Group_Pool_branch_size_mean ML_Pool_branch_size_mean];
%     Group_Pool_thresholded_branch_number_mean=[Group_Pool_thresholded_branch_number_mean ML_Pool_thresholded_branch_number_mean];
%     Group_Pool_branch_cellmovement_std=[Group_Pool_branch_cellmovement_std ML_Pool_branch_cellmovement_std];
%     Group_Pool_fila_branch_orientation_pool_std=[Group_Pool_fila_branch_orientation_pool_std ML_Pool_fila_branch_orientation_pool_std];
%     Group_Pool_fila_trajectory_orientation_pool_std=[Group_Pool_fila_trajectory_orientation_pool_std ML_Pool_fila_trajectory_orientation_pool_std];
%     Group_Pool_fila_branch_orientation_pool=[Group_Pool_fila_branch_orientation_pool; ML_Pool_fila_branch_orientation_pool];
%     Group_Pool_fila_trajectory_orientation_pool=[Group_Pool_fila_trajectory_orientation_pool; ML_Pool_fila_trajectory_orientation_pool];
%     Group_Pool_fila_branch_orientation_pool_slow=[Group_Pool_fila_branch_orientation_pool_slow; ML_Pool_fila_branch_orientation_pool_slow];
%     Group_Pool_fila_trajectory_orientation_pool_slow=[Group_Pool_fila_trajectory_orientation_pool_slow; ML_Pool_fila_trajectory_orientation_pool_slow];
%     Group_Pool_fila_branch_orientation_pool_fast=[Group_Pool_fila_branch_orientation_pool_fast; ML_Pool_fila_branch_orientation_pool_fast];
%     Group_Pool_fila_trajectory_orientation_pool_fast=[Group_Pool_fila_trajectory_orientation_pool_fast; ML_Pool_fila_trajectory_orientation_pool_fast];
%     Group_Pool_CompletedFrame_last =[Group_Pool_CompletedFrame_last ML_Pool_CompletedFrame_last];
%     
%     
%     Group_Pool_whole_cell_vim_seg_total = [Group_Pool_whole_cell_vim_seg_total; ML_Pool_whole_cell_vim_seg_total];
%     Group_Pool_whole_cell_vim_seg_mean  = [Group_Pool_whole_cell_vim_seg_mean;  ML_Pool_whole_cell_vim_seg_mean ];
%     Group_Pool_whole_cell_vim_nms_total = [Group_Pool_whole_cell_vim_nms_total; ML_Pool_whole_cell_vim_nms_total];
%     Group_Pool_whole_cell_vim_nms_mean  = [Group_Pool_whole_cell_vim_nms_mean; ML_Pool_whole_cell_vim_nms_mean ];
%     
%     Group_Pool_branch_seg_total = [Group_Pool_branch_seg_total ML_Pool_branch_seg_total];
%     Group_Pool_branch_seg_mean  = [Group_Pool_branch_seg_mean  ML_Pool_branch_seg_mean ];
%     Group_Pool_branch_nms_total = [Group_Pool_branch_nms_total ML_Pool_branch_nms_total];
%     Group_Pool_branch_nms_mean  = [Group_Pool_branch_nms_mean  ML_Pool_branch_nms_mean ];
%    
%     for iShorttime = 1 : numel(short_time_frames)        
%         Group_Pool_Travel_Speed_short_time{iShorttime} = [ Group_Pool_Travel_Speed_short_time{iShorttime};   ML_Pool_Travel_Speed_short_time{iShorttime}];
%         Group_Pool_vim_fila_shorttime{iShorttime} = [ Group_Pool_vim_fila_shorttime{iShorttime};    ML_Pool_vim_fila_shorttime{iShorttime};];
%         Group_Pool_vim_nms_shorttime{iShorttime} = [ Group_Pool_vim_nms_shorttime{iShorttime};    ML_Pool_vim_nms_shorttime{iShorttime};];
%         Group_Pool_vim_int_shorttime{iShorttime} = [ Group_Pool_vim_int_shorttime{iShorttime};    ML_Pool_vim_nms_shorttime{iShorttime};];        
%         Group_Pool_persistency_shorttime{iShorttime} = [ Group_Pool_persistency_shorttime{iShorttime};    ML_Pool_persistency_shorttime{iShorttime};];        
%     end
%         
%     for iB = 1 : numel(branch_percent_threshold)
%         Group_Pool_branch_number_weighted{iB} = [ Group_Pool_branch_number_weighted{iB}; ML_Pool_branch_number_weighted{iB};];
%     end
%     
%     
%     for iSpeed = 1 : numel(speed_percent_threshold)
%         for iPvalue = 1 : numel(branch_orientation_p_threshold)
%             for iShorttime = 1 : numel(short_time_frames)                
%                 Group_Pool_reorienting_status_short_time{iShorttime,iSpeed,iPvalue} = ...
%                     [Group_Pool_reorienting_status_short_time{iShorttime,iSpeed,iPvalue}; ...
%                     ML_Pool_reorienting_status_short_time{iShorttime,iSpeed,iPvalue};];                
%             end
%         end
%     end

      Group_Pool_curvature_map_pool= [Group_Pool_curvature_map_pool;ML_Pool_curvature_map_pool];
    Group_Pool_filament_alignment_map_pool= [Group_Pool_filament_alignment_map_pool;ML_Pool_filament_alignment_map_pool];
    Group_Pool_filament_density_pool= [Group_Pool_filament_density_pool;ML_Pool_filament_density_pool];
    Group_Pool_prot_or_retr_fila_pool= [Group_Pool_prot_or_retr_fila_pool;ML_Pool_prot_or_retr_fila_pool];
    Group_Pool_curvature_map_full_pool= [Group_Pool_curvature_map_full_pool;ML_Pool_curvature_map_full_pool];
    Group_Pool_filament_alignment_map_full_pool= [Group_Pool_filament_alignment_map_full_pool;ML_Pool_filament_alignment_map_full_pool];
    
    Group_Pool_prot_or_retr_fila_full_pool= [Group_Pool_prot_or_retr_fila_full_pool;ML_Pool_prot_or_retr_fila_full_pool];
    Group_Pool_scale_map_pool= [Group_Pool_scale_map_pool;ML_Pool_scale_map_pool];
    
    Group_Pool_last_curvature_map_pool= [Group_Pool_last_curvature_map_pool;ML_Pool_last_curvature_map_pool];
    Group_Pool_last_filament_alignment_map_pool= [Group_Pool_last_filament_alignment_map_pool;ML_Pool_last_filament_alignment_map_pool];
    Group_Pool_last_filament_density_pool= [Group_Pool_last_filament_density_pool;ML_Pool_last_filament_density_pool];
    Group_Pool_last_prot_or_retr_fila_pool= [Group_Pool_last_prot_or_retr_fila_pool;ML_Pool_last_prot_or_retr_fila_pool];
    Group_Pool_last_curvature_map_full_pool= [Group_Pool_last_curvature_map_full_pool;ML_Pool_last_curvature_map_full_pool];
    Group_Pool_last_filament_alignment_map_full_pool= [Group_Pool_last_filament_alignment_map_full_pool;ML_Pool_last_filament_alignment_map_full_pool];
    Group_Pool_last_filament_density_full_pool= [Group_Pool_last_filament_density_full_pool;ML_Pool_last_filament_density_full_pool];
    Group_Pool_last_scale_map_pool= [Group_Pool_last_scale_map_pool;ML_Pool_last_scale_map_pool];
    
        
    Group_Pool_int_map_pool= [Group_Pool_int_map_pool;ML_Pool_int_map_pool];
    Group_Pool_st_map_pool= [Group_Pool_st_map_pool;ML_Pool_st_map_pool];
    Group_Pool_int_map_full_pool= [Group_Pool_int_map_full_pool;ML_Pool_int_map_full_pool];
    Group_Pool_st_map_full_pool= [Group_Pool_st_map_full_pool;ML_Pool_st_map_full_pool];
    
    Group_Pool_last_int_map_pool= [Group_Pool_last_int_map_pool;ML_Pool_last_int_map_pool];
    Group_Pool_last_st_map_pool= [Group_Pool_last_st_map_pool;ML_Pool_last_st_map_pool];
    Group_Pool_last_int_map_full_pool= [Group_Pool_last_int_map_full_pool;ML_Pool_last_int_map_full_pool];
    Group_Pool_last_st_map_full_pool= [Group_Pool_last_st_map_full_pool;ML_Pool_last_st_map_full_pool];
        
    
    save([ML_name,'Group_prot_retract_Pools.mat']); 
end

%% Plotting for all MLs
colorarray = rand(Group_Pool_last_st_map_full_pool,3);

  h16 = figure(16);
        subplot(3,2,1);
        hist(BA_output.curvature_map_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('Curvature full prot-ed');
        subplot(3,2,2);
        hist(BA_output.curvature_map_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('Curvature full retract-ed');
        subplot(3,2,3);
        hist(BA_output.filament_alignment_map_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('Filament Alignment full prot-ed');
        subplot(3,2,4);
        hist(BA_output.filament_alignment_map_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('Filament Alignment full retract-ed');
        subplot(3,2,5);
        hist(BA_output.filament_density_full_pool(BA_output.prot_or_retr_fila_full_pool==1),30);
        title('Filament density full prot-ed');
        subplot(3,2,6);
        hist(BA_output.filament_density_full_pool(BA_output.prot_or_retr_fila_full_pool==-1),30);
        title('Filament density full retract-ed');
        saveas(h16,[outputPath,filesep,'this_filament_prot_retract_full.tif']);
        saveas(h16,[outputPath,filesep,'this_filament_prot_retract_full.fig']);
        
%% % plot the this group result
h3 = figure(3);
plot(Group_Pool_branch_duration_array, Group_Pool_branch_vif_mean_intensity,'.');
xlabel(['Branch Duration,sample size:',num2str(numel(Group_Pool_branch_duration_array))],'Fontsize',13);
ylabel('Branch Vim mean Int','Fontsize',13);
title('Branch: Duration vs Vim Mean Int','Fontsize',13);
set(gca,'fontsize',13)
saveas(h3,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.fig']);
saveas(h3,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.tif']);
print(h3,'-depsc',[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.eps']);

