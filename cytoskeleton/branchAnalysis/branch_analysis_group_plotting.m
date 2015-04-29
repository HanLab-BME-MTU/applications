function BA_Group = branch_analysis_group_plotting(ML_name_cell, T_branchsize, T_branchduration, Group_ROOT_DIR, branch_percent_threshold, short_time_frames, speed_percent_threshold,branch_orientation_p_threshold)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML_array:     The array of movieList objects in this group to be analyzed

nList = length(ML_name_cell);

BA_Group = cell(1, nList);

if(nargin<2)
    T_branchsize=200;
end
if(nargin<3)
    T_branchduration=2;
end
if(nargin<4)
    Group_ROOT_DIR=[];
end
if(nargin<5)
    branch_percent_threshold = 1; % anything larger than 1% of cell size
end
if(nargin<6)
    short_time_frames = 1; % check things on each period of time  1 frame
end
if(nargin<7)
    speed_percent_threshold = 1; % the speed to be considered as valid, not reorienting, 1% of cell demension( sqrt of cell pixels) is very low
end
if(nargin<8)
    branch_orientation_p_threshold = 0.999; %  the threshold to consider branch orientation as different 0.99 means they can be really similar to each other
end

T_speed=4;

close all;

for iML = 1 : nList
    
    ML_name = ML_name_cell{iML};
    
    load(ML_name);
    
    ML_ROOT_DIR = ML.outputDirectory_;
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    % if the batch results exist, load it
    if 0
%         (exist([ML_ROOT_DIR,filesep,'movieList_BA_results_gathered_balloon.mat'], 'file'))
%         load([ML_ROOT_DIR,filesep,'movieList_BA_results_gathered_balloon.mat'],'BA_output_ML_cell');
    else
        % otherwise load one by one
        BA_output_ML_cell= cell(1,1,1);
        
        for iM  = 1 :movieNumber
            try
                % load this movie
                load(ML.movieDataFile_{iM});
            catch
                disp('The MD file is missing');
                continue;
            end
            % Now MD in workspace
            ROOT_DIR = MD.outputDirectory_;
            
            % the number of channels
            nChannel = length(MD.channels_);
            
            for iChannel = 1 : nChannel
                for iCell = 1 : 20
                    display(['Checking: iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                    
                    % the folder name if there is marking
                    outputPath = [ROOT_DIR,filesep,'BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];
                    
                    % check if this folder exist
                    if(exist(outputPath,'dir'))
                        % if it exist, try to do the branch analysis
                        try
                            if(exist([outputPath,filesep,'branch_analysis_results_balloon.mat'],'file'))
                                load([outputPath,filesep,'branch_analysis_results_balloon.mat'],'BA_output');
                            else
                                continue;
                            end
                            
                            display(['Found: iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                            Identifier = ['ML',num2str(iML),'M',num2str(iM),'C',num2str(iCell)];
                            BA_output.Identifier = Identifier;
                            % save
                            BA_output_ML_cell{1, iM}{iChannel, iCell} = BA_output;
                        end
                    end
                end
            end
        end
    end
    
    BA_Group{1, iML} = BA_output_ML_cell;
end



% initialize the pools

Group_Pool_Travel_Length = [];
Group_Pool_Travel_Distance = [];
Group_Pool_Travel_Speed = [];
Group_Pool_Travel_Speed_without_pausing = [];
Group_Pool_Cell_Marked_Frame_Number = [];
Group_Pool_branch_number_tracked = [];
Group_Pool_branch_number_max = [];
Group_Pool_branch_number_mean = [];
Group_Pool_branch_duration_array = [];
Group_Pool_branch_duration_mean= [];
Group_Pool_branch_vif_mean_intensity= [];
Group_Pool_protrusion_vif_mean_intensity= [];
Group_Pool_retraction_vif_mean_intensity= [];
Group_Pool_whole_cell_vif_mean_intensity= [];
Group_Pool_whole_cell_size_mean= [];
Group_Pool_branch_number_mean_pat=[];
Group_Pool_whole_cell_vim_totalamount_mean=[];
Group_Pool_whole_cell_vif_mean_intensity_pat=[];
Group_Pool_branch_size_mean=[];
Group_Pool_thresholded_branch_number_mean=[];
Group_Pool_branch_cellmovement_std=[];
Group_Pool_fila_branch_orientation_pool_std=[];
Group_Pool_fila_trajectory_orientation_pool_std=[];
Group_Pool_fila_branch_orientation_pool=[];
Group_Pool_fila_trajectory_orientation_pool=[];
Group_Pool_fila_branch_orientation_pool_slow=[];
Group_Pool_fila_trajectory_orientation_pool_slow=[];
Group_Pool_fila_branch_orientation_pool_fast=[];
Group_Pool_fila_trajectory_orientation_pool_fast=[];
Group_Pool_pool_vif_intensity=[];
Group_Pool_CompletedFrame_last=[];
Group_Pool_whole_cell_vim_seg_total = [];
Group_Pool_whole_cell_vim_seg_mean  = [];
Group_Pool_whole_cell_vim_nms_total = [];
Group_Pool_whole_cell_vim_nms_mean  = [];
Group_Pool_branch_seg_total = [];
Group_Pool_branch_seg_mean  = [];
Group_Pool_branch_nms_total = [];
Group_Pool_branch_nms_mean  = [];
Group_Pool_branch_number_weighted = cell(1,numel(branch_percent_threshold));
Group_Pool_Travel_Speed_short_time = cell(1, numel(short_time_frames));
Group_Pool_reorienting_status_short_time = cell{numel(short_time_frames),numel(speed_percent_threshold),numel(branch_orientation_p_threshold));
Group_Pool_vim_fila_shorttime = cell(1, numel(short_time_frames));
Group_Pool_vim_nms_shorttime = cell(1, numel(short_time_frames));
Group_Pool_vim_int_shorttime = cell(1, numel(short_time_frames));

Identifier_cell = [];

iAllCell = 0;

% if no input, the DIR of last movieList
if(isempty(Group_ROOT_DIR))
    Group_ROOT_DIR = ML_ROOT_DIR;
end

if(exist(Group_ROOT_DIR,'dir')==0)
    mkdir(Group_ROOT_DIR);
end

%% For each ML
for iML = 1 : nList
    BA_output_ML_cell = BA_Group{1, iML};
    
    ML_Pool_Travel_Length = [];
    ML_Pool_Travel_Distance = [];
    ML_Pool_Travel_Speed = [];
    ML_Pool_Cell_Marked_Frame_Number = [];
    ML_Pool_branch_number_tracked = [];
    ML_Pool_branch_number_max = [];
    ML_Pool_branch_number_mean = [];
    ML_Pool_branch_duration_array = [];
    ML_Pool_branch_duration_mean= [];
    ML_Pool_branch_vif_mean_intensity= [];
    ML_Pool_protrusion_vif_mean_intensity= [];
    ML_Pool_retraction_vif_mean_intensity= [];
    ML_Pool_whole_cell_vif_mean_intensity= [];
    ML_Pool_whole_cell_size_mean= [];
    ML_Pool_branch_number_mean_pat=[];
    ML_Pool_whole_cell_vim_totalamount_mean=[];
    ML_Pool_whole_cell_vif_mean_intensity_pat=[];
    ML_Pool_branch_size_mean=[];
    ML_Pool_thresholded_branch_number_mean=[];
    ML_Pool_branch_cellmovement_std=[];
    ML_Pool_fila_branch_orientation_pool_std=[];
    ML_Pool_fila_trajectory_orientation_pool_std=[];
    ML_Pool_fila_branch_orientation_pool=[];
    ML_Pool_fila_trajectory_orientation_pool=[];
    ML_Pool_fila_branch_orientation_pool_slow=[];
    ML_Pool_fila_trajectory_orientation_pool_slow=[];
    ML_Pool_fila_branch_orientation_pool_fast=[];
    ML_Pool_fila_trajectory_orientation_pool_fast=[];
    ML_Pool_CompletedFrame_last = [];
    ML_Pool_pool_vif_intensity=[];
    ML_Pool_whole_cell_vim_seg_total = [];
    ML_Pool_whole_cell_vim_seg_mean = [];
    ML_Pool_whole_cell_vim_nms_total = [];
    ML_Pool_whole_cell_vim_nms_mean = [];
    ML_Pool_branch_seg_total = [];
    ML_Pool_branch_seg_mean = [];
    ML_Pool_branch_nms_total = [];
    ML_Pool_branch_nms_mean = [];
    ML_Pool_Travel_Speed_without_pausing = [];
    
    ML_Pool_branch_number_weighted = cell(1,numel(branch_percent_threshold));
    ML_Pool_Travel_Speed_short_time = cell{1, numel(short_time_frames)};
    ML_Pool_reorienting_status_short_time = cell{numel(short_time_frames),numel(speed_percent_threshold),numel(branch_orientation_p_threshold)};
    ML_Pool_vim_fila_shorttime = cell{1, numel(short_time_frames)};
    ML_Pool_vim_nms_shorttime = cell{1, numel(short_time_frames)};
    ML_Pool_vim_int_shorttime = cell{1, numel(short_time_frames)};
    
    for iB = 1 : numel(branch_percent_threshold)
        
        percent_threshold = branch_percent_threshold(iB);
        
        branch_size_matrix = BA_output.branch_size_matrix;
        cell_size_pool = BA_output.cell_size_array;
        
        cell_size_matrix = repmat(cell_size_pool, 1,size(branch_size_matrix,2));
        
        
        branch_size_percent_matrix = branch_size_matrix./cell_size_matrix*100;
        
        branch_size_percent_thresholded_matrix  = branch_size_percent_matrix/percent_threshold;
        branch_size_percent_thresholded_matrix(branch_size_percent_thresholded_matrix>1)=1;
        
        ML_Pool_branch_number_weighted{iB} = [ ML_Pool_branch_number_weighted{iB}; nanmean(branch_size_percent_thresholded_matrix,1)];
        
    end
       
   movieNumber =  length(ML.movieDataFile_);
        
    for iM  = 1 : movieNumber
        nChannel = 2;
        for iChannel = 1 : nChannel
            for iCell = 1 :  10
                try
                    BA_output = BA_output_ML_cell{1, iM}{iChannel, iCell};
                catch
                    continue;
                end
                
                if(~isempty(BA_output))
                    ML_Pool_CompletedFrame_last=[ML_Pool_CompletedFrame_last BA_output.CompletedFrame(end)];    
                    ML_Pool_branch_size_mean = [ML_Pool_branch_size_mean BA_output.branch_mean_size];
                    ML_Pool_whole_cell_size_mean= [ML_Pool_whole_cell_size_mean BA_output.whole_cell_size_mean];
                    ML_Pool_branch_number_mean_pat = [ML_Pool_branch_number_mean_pat repmat(BA_output.branch_number_mean, [1, length(BA_output.branch_vif_mean_intensity)])];
                    Identifier_cell{1,numel(Identifier_cell)+1} = BA_output.Identifier;
                    
                    Thresholded_branch_number_mean = sum(BA_output.branch_duration_array(find(BA_output.branch_mean_size>T_branchsize & BA_output.branch_duration_array>T_branchduration)))...
                        ./ BA_output.cell_marked_frame_number;
                    ML_Pool_thresholded_branch_number_mean = [ML_Pool_thresholded_branch_number_mean Thresholded_branch_number_mean];
                    
                    
                    ML_Pool_Travel_Length = [ML_Pool_Travel_Length BA_output.cell_travel_length];
                    ML_Pool_Travel_Distance = [ML_Pool_Travel_Distance BA_output.cell_travel_distance];
                    ML_Pool_Travel_Speed = [ML_Pool_Travel_Speed BA_output.cell_travel_speed];
                    ML_Pool_Cell_Marked_Frame_Number = [ML_Pool_Cell_Marked_Frame_Number BA_output.cell_marked_frame_number];
                    ML_Pool_branch_number_tracked = [ML_Pool_branch_number_tracked BA_output.branch_number_tracked];
                    ML_Pool_branch_number_max = [ML_Pool_branch_number_max BA_output.branch_number_max];
                    ML_Pool_branch_number_mean = [ML_Pool_branch_number_mean BA_output.branch_number_mean];
                    ML_Pool_branch_duration_array = [ML_Pool_branch_duration_array  BA_output.branch_duration_array];
                    ML_Pool_branch_duration_mean= [ML_Pool_branch_duration_mean  BA_output.branch_duration_mean];
                    
                    ML_Pool_whole_cell_vif_mean_intensity_pat = [ML_Pool_whole_cell_vif_mean_intensity_pat repmat(BA_output.whole_cell_vif_mean_intensity, [1, length(BA_output.branch_vif_mean_intensity)])];
                    ML_Pool_whole_cell_vim_totalamount_mean = [ML_Pool_whole_cell_vim_totalamount_mean BA_output.whole_cell_vim_totalamount_mean];
                    ML_Pool_branch_vif_mean_intensity= [ML_Pool_branch_vif_mean_intensity BA_output.branch_vif_mean_intensity];
                    ML_Pool_protrusion_vif_mean_intensity= [ML_Pool_protrusion_vif_mean_intensity BA_output.protrusion_vif_mean_intensity];
                    ML_Pool_retraction_vif_mean_intensity= [ML_Pool_retraction_vif_mean_intensity BA_output.retraction_vif_mean_intensity];
                    ML_Pool_whole_cell_vif_mean_intensity= [ML_Pool_whole_cell_vif_mean_intensity BA_output.whole_cell_vif_mean_intensity];
                    
                    ML_Pool_pool_vif_intensity= [ML_Pool_pool_vif_intensity; BA_output.pool_all_vif_intensity];
                    
                    
                    ML_Pool_whole_cell_vim_seg_total = [ML_Pool_whole_cell_vim_seg_total; BA_output.whole_cell_vim_seg_total];
                    ML_Pool_whole_cell_vim_seg_mean = [ML_Pool_whole_cell_vim_seg_mean; BA_output.whole_cell_vim_seg_mean];
                    ML_Pool_whole_cell_vim_nms_total = [ML_Pool_whole_cell_vim_nms_total; BA_output.whole_cell_vim_nms_total];
                    ML_Pool_whole_cell_vim_nms_mean = [ML_Pool_whole_cell_vim_nms_mean; BA_output.whole_cell_vim_nms_mean];
                    
                    ML_Pool_branch_seg_total = [ML_Pool_branch_seg_total BA_output.branch_seg_total];
                    ML_Pool_branch_seg_mean = [ML_Pool_branch_seg_mean BA_output.branch_seg_mean];
                    ML_Pool_branch_nms_total = [ML_Pool_branch_nms_total BA_output.branch_nms_total];
                    ML_Pool_branch_nms_mean = [ML_Pool_branch_nms_mean BA_output.branch_nms_mean];
                    
                    if isfield(BA_output,'speed_marked_frames')
                        S = BA_output.speed_marked_frames;
                        ML_Pool_Travel_Speed_without_pausing = [ML_Pool_Travel_Speed_without_pausing mean(S(S>T_speed))];
                        smoothed_speed_marked_frames = BA_output.smoothed_speed_marked_frames;
                        center_x = BA_output.center_x;
                        center_y = BA_output.center_y;
                        center_x = center_x(:);
                        center_y = center_y(:);
                        % assume the last frame is moving in the same speed
                        % as the end-1 to end
                        smoothed_speed_marked_frames(numel(smooth_center_x)) = smoothed_speed_marked_frames(numel(smooth_center_x)-1);
                        smoothed_speed_marked_frames = smoothed_speed_marked_frames(:);
      
      
                        %%  speed with shorttime definition
                        for iShorttime = 1 : numel(short_time_frames)
                            ShortTimeFrames = short_time_frames(iShorttime);
%                             give_up_frame = floor(ShortTimeFrames/2);
                            give_up_frame = 1;
                            H = ones(ShortTimeFrames,1)/ShortTimeFrames;
                            
                            smoothed_speed_marked_frames_shorttime = imfilter(smoothed_speed_marked_frames,H,'same','replicate');
                            
                            vif_fila_mean = BA_output.cell_vif_seg_total_array ./(BA_output.cell_size_array);
                            vif_nms_mean = BA_output.cell_vif_nms_total_array ./(BA_output.cell_size_array);
                            vif_int_mean = BA_output.cell_vimtotal_pool ./(BA_output.cell_size_array);
                            
                            vif_fila_mean = vif_fila_mean(:);
                            vif_nms_mean  = vif_nms_mean(:);
                            vif_int_mean  = vif_int_mean(:);
                            
                            vif_fila_mean_shorttime = imfilter(vif_fila_mean,H,'same','replicate');
                            vif_nms_mean_shorttime = imfilter(vif_nms_mean,H,'same','replicate');
                            vif_int_mean_shorttime = imfilter(vif_int_mean,H,'same','replicate');
                            
                            ML_Pool_Travel_Speed_short_time{iShorttime} = [ML_Pool_Travel_Speed_short_time{iShorttime}; smoothed_speed_marked_frames_shorttime(give_up_frame:end-give_up_frame+1,1)] ;
                            ML_Pool_vim_fila_shorttime{iShorttime} = [ML_Pool_vim_fila_shorttime{iShorttime}; vif_fila_mean_shorttime(give_up_frame:end-give_up_frame+1,1)];
                            ML_Pool_vim_nms_shorttime{iShorttime} = [ML_Pool_vim_nms_shorttime{iShorttime}; vif_nms_mean_shorttime(give_up_frame:end-give_up_frame+1,1)];
                            ML_Pool_vim_int_shorttime{iShorttime} = [ML_Pool_vim_int_shorttime{iShorttime}; vif_int_mean_shorttime(give_up_frame:end-give_up_frame+1,1)];                            
                        end
                        
                        %% branch number with weight
                        for iB = 1 : numel(branch_percent_threshold)
                            branch_size_percent = BA_output.branch_size_matrix./(repmat(BA_output.cell_size_array, 1, size(BA_output.branch_size_matrix,2)))*100;
                            branch_size_percent(branch_size_percent> branch_percent_threshold(iB)) = branch_percent_threshold(iB);
                            branch_size_percent = branch_size_percent/branch_size_percent(iB);
                            branch_number_afterweight = sum(branch_size_percent,2);
                            ML_Pool_branch_number_weighted{iB} = [ ML_Pool_branch_number_weighted{iB}; branch_number_afterweight];
                        end
                        %%
    
      for iSpeed = 1 : numel(speed_percent_threshold)
          for iPvalue = 1 : numel(speed_percent_threshold)
              for iShorttime = 1 : numel(short_time_frames)
                  ShortTimeFrames = short_time_frames(iShorttime);
                 
%%
        
              end
          end
      end
      
      
      
      
   end
                    
                    
                    % some old analysis doens't have orient std,
                    % but as space holder, add nan
                    try
                        ML_Pool_branch_cellmovement_std= [ML_Pool_branch_cellmovement_std BA_output.branch_cellmovement_std];
                        
                    catch
                        ML_Pool_branch_cellmovement_std= [ML_Pool_branch_cellmovement_std nan];
                    end
                    
                    % some old analysis doens't have orient std,
                    % but as space holder, add nan
                    try
                        ML_Pool_fila_branch_orientation_pool_std= [ML_Pool_fila_branch_orientation_pool_std BA_output.fila_branch_orientation_pool_std];
                        
                    catch
                        ML_Pool_fila_branch_orientation_pool_std= [ML_Pool_fila_branch_orientation_pool_std nan];
                    end
                    
                    try
                        ML_Pool_fila_trajectory_orientation_pool_std= [ML_Pool_fila_trajectory_orientation_pool_std BA_output.fila_trajectory_orientation_pool_std];
                        
                    catch
                        ML_Pool_fila_trajectory_orientation_pool_std= [ML_Pool_fila_trajectory_orientation_pool_std nan];
                    end
                    
                    %some old analysis doesn't have these pools
                    try
                        ML_Pool_branch_cellmovement_pool= [ML_Pool_branch_cellmovement_pool; BA_output.branch_trajectory_orientation_pool];
                        
                        if(BA_output.cell_travel_speed>7)
                            ML_Pool_branch_cellmovement_pool_fast= [ML_Pool_branch_cellmovement_pool_fast; BA_output.branch_trajectory_orientation_pool];
                        else
                            ML_Pool_branch_cellmovement_pool_slow= [ML_Pool_branch_cellmovement_pool_slow; BA_output.branch_trajectory_orientation_pool];
                        end
                    end
                    
                    
                    
                    try
                        ML_Pool_fila_branch_orientation_pool= [ML_Pool_fila_branch_orientation_pool; BA_output.fila_branch_orientation_pool];
                        
                        if(BA_output.cell_travel_speed>7)
                            ML_Pool_fila_branch_orientation_pool_fast=[ML_Pool_fila_branch_orientation_pool_fast; BA_output.fila_branch_orientation_pool];
                        else
                            ML_Pool_fila_branch_orientation_pool_slow=[ML_Pool_fila_trajectory_orientation_pool_slow; BA_output.fila_trajectory_orientation_pool];
                        end
                    end
                    
                    
                    try
                        ML_Pool_fila_trajectory_orientation_pool= [ML_Pool_fila_trajectory_orientation_pool; BA_output.fila_trajectory_orientation_pool];
                        
                        if(BA_output.cell_travel_speed>7)
                            ML_Pool_fila_trajectory_orientation_pool_fast=[ML_Pool_fila_trajectory_orientation_pool_fast; BA_output.fila_trajectory_orientation_pool];
                        else
                            ML_Pool_fila_trajectory_orientation_pool_slow=[ML_Pool_fila_branch_orientation_pool_slow; BA_output.fila_branch_orientation_pool];
                        end
                    end
                    iAllCell = iAllCell + 1;
                end
            end
        end
    end
    
    %% vim normalization
    vim_max = prctile(double(ML_Pool_pool_vif_intensity),98);
    
    display('50')
    prctile(double(ML_Pool_pool_vif_intensity),50)
    display('75')
    prctile(double(ML_Pool_pool_vif_intensity),75)
    display('95')
    prctile(double(ML_Pool_pool_vif_intensity),95)
    display('98')
    prctile(double(ML_Pool_pool_vif_intensity),98)
    display('99')
    prctile(double(ML_Pool_pool_vif_intensity),99)
    
    
    ML_Pool_branch_vif_mean_intensity = ML_Pool_branch_vif_mean_intensity/vim_max;
    ML_Pool_protrusion_vif_mean_intensity = ML_Pool_protrusion_vif_mean_intensity/vim_max;
    ML_Pool_retraction_vif_mean_intensity = ML_Pool_retraction_vif_mean_intensity/vim_max;
    ML_Pool_whole_cell_vif_mean_intensity = ML_Pool_whole_cell_vif_mean_intensity/vim_max;
    ML_Pool_whole_cell_vif_mean_intensity_pat = ML_Pool_whole_cell_vif_mean_intensity_pat/vim_max;
    ML_Pool_whole_cell_vim_totalamount_mean = ML_Pool_whole_cell_vim_totalamount_mean/vim_max;
    
    
    for shorttime also need this normalization
        """""""""""""""""""""""""""""""
   
    
    %% in the end put data from this ML to the pool of all MLs
        
    Group_Pool_Travel_Length = [Group_Pool_Travel_Length ML_Pool_Travel_Length];
    Group_Pool_Travel_Distance = [Group_Pool_Travel_Distance ML_Pool_Travel_Distance];
    Group_Pool_Travel_Speed = [Group_Pool_Travel_Speed ML_Pool_Travel_Speed];
    Group_Pool_Travel_Speed_without_pausing = [Group_Pool_Travel_Speed_without_pausing ML_Pool_Travel_Speed_without_pausing];
    Group_Pool_Cell_Marked_Frame_Number = [Group_Pool_Cell_Marked_Frame_Number ML_Pool_Cell_Marked_Frame_Number];
    Group_Pool_branch_number_tracked = [Group_Pool_branch_number_tracked ML_Pool_branch_number_tracked];
    Group_Pool_branch_number_max = [Group_Pool_branch_number_max ML_Pool_branch_number_max];
    Group_Pool_branch_number_mean = [Group_Pool_branch_number_mean ML_Pool_branch_number_mean];
    Group_Pool_branch_duration_array = [Group_Pool_branch_duration_array ML_Pool_branch_duration_array];
    Group_Pool_branch_duration_mean= [Group_Pool_branch_duration_mean ML_Pool_branch_duration_mean];
    Group_Pool_branch_vif_mean_intensity= [Group_Pool_branch_vif_mean_intensity ML_Pool_branch_vif_mean_intensity];
    Group_Pool_protrusion_vif_mean_intensity= [Group_Pool_protrusion_vif_mean_intensity ML_Pool_protrusion_vif_mean_intensity];
    Group_Pool_retraction_vif_mean_intensity= [Group_Pool_retraction_vif_mean_intensity ML_Pool_retraction_vif_mean_intensity];
    Group_Pool_whole_cell_vif_mean_intensity= [Group_Pool_whole_cell_vif_mean_intensity ML_Pool_whole_cell_vif_mean_intensity];
    Group_Pool_whole_cell_size_mean= [Group_Pool_whole_cell_size_mean ML_Pool_whole_cell_size_mean];
    Group_Pool_branch_number_mean_pat=[Group_Pool_branch_number_mean_pat ML_Pool_branch_number_mean_pat];
    Group_Pool_whole_cell_vim_totalamount_mean=[Group_Pool_whole_cell_vim_totalamount_mean ML_Pool_whole_cell_vim_totalamount_mean];
    Group_Pool_whole_cell_vif_mean_intensity_pat=[Group_Pool_whole_cell_vif_mean_intensity_pat ML_Pool_whole_cell_vif_mean_intensity_pat];
    Group_Pool_branch_size_mean=[Group_Pool_branch_size_mean ML_Pool_branch_size_mean];
    Group_Pool_thresholded_branch_number_mean=[Group_Pool_thresholded_branch_number_mean ML_Pool_thresholded_branch_number_mean];
    Group_Pool_branch_cellmovement_std=[Group_Pool_branch_cellmovement_std ML_Pool_branch_cellmovement_std];
    Group_Pool_fila_branch_orientation_pool_std=[Group_Pool_fila_branch_orientation_pool_std ML_Pool_fila_branch_orientation_pool_std];
    Group_Pool_fila_trajectory_orientation_pool_std=[Group_Pool_fila_trajectory_orientation_pool_std ML_Pool_fila_trajectory_orientation_pool_std];
    Group_Pool_fila_branch_orientation_pool=[Group_Pool_fila_branch_orientation_pool; ML_Pool_fila_branch_orientation_pool];
    Group_Pool_fila_trajectory_orientation_pool=[Group_Pool_fila_trajectory_orientation_pool; ML_Pool_fila_trajectory_orientation_pool];
    Group_Pool_fila_branch_orientation_pool_slow=[Group_Pool_fila_branch_orientation_pool_slow; ML_Pool_fila_branch_orientation_pool_slow];
    Group_Pool_fila_trajectory_orientation_pool_slow=[Group_Pool_fila_trajectory_orientation_pool_slow; ML_Pool_fila_trajectory_orientation_pool_slow];
    Group_Pool_fila_branch_orientation_pool_fast=[Group_Pool_fila_branch_orientation_pool_fast; ML_Pool_fila_branch_orientation_pool_fast];
    Group_Pool_fila_trajectory_orientation_pool_fast=[Group_Pool_fila_trajectory_orientation_pool_fast; ML_Pool_fila_trajectory_orientation_pool_fast];
    Group_Pool_CompletedFrame_last =[Group_Pool_CompletedFrame_last ML_Pool_CompletedFrame_last];
    
    
    Group_Pool_whole_cell_vim_seg_total = [Group_Pool_whole_cell_vim_seg_total; ML_Pool_whole_cell_vim_seg_total];
    Group_Pool_whole_cell_vim_seg_mean  = [Group_Pool_whole_cell_vim_seg_mean;  ML_Pool_whole_cell_vim_seg_mean ];
    Group_Pool_whole_cell_vim_nms_total = [Group_Pool_whole_cell_vim_nms_total; ML_Pool_whole_cell_vim_nms_total];
    Group_Pool_whole_cell_vim_nms_mean  = [Group_Pool_whole_cell_vim_nms_mean; ML_Pool_whole_cell_vim_nms_mean ];
    
    Group_Pool_branch_seg_total = [Group_Pool_branch_seg_total ML_Pool_branch_seg_total];
    Group_Pool_branch_seg_mean  = [Group_Pool_branch_seg_mean  ML_Pool_branch_seg_mean ];
    Group_Pool_branch_nms_total = [Group_Pool_branch_nms_total ML_Pool_branch_nms_total];
    Group_Pool_branch_nms_mean  = [Group_Pool_branch_nms_mean  ML_Pool_branch_nms_mean ];
    
end

%% Plotting for all MLs
Group_Pool_tracked_branch_d_frames = Group_Pool_branch_number_tracked./Group_Pool_Cell_Marked_Frame_Number;

colorarray = rand(numel(Group_Pool_branch_vif_mean_intensity),3);



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


size(Group_Pool_protrusion_vif_mean_intensity)

h4 = figure(4);
[h_p,bin] = hist(Group_Pool_protrusion_vif_mean_intensity,0.0:0.03:1.1);
[h_r,bin] = hist(Group_Pool_retraction_vif_mean_intensity,0.0:0.03:1.1);
bar([reshape(bin, 1, numel(bin)); reshape(bin, 1, numel(bin))]',[reshape(h_p, 1, numel(h_p));reshape(h_r, 1, numel(h_r))]');
legend('Protrusion Region Vim','Retraction Region Vim');
axis([0.00 0.30 0 max([h_p h_r])])
set(gca,'fontsize',13)
saveas(h4,[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.fig']);
saveas(h4,[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.tif']);
print(h4,'-depsc',[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.eps']);

% mean vim vs speed
h7 = figure(7); hold off;
%  plot(Group_Pool_branch_number_mean, Group_Pool_Travel_Speed,'o');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_Travel_Speed(i)+0.2, num2str(i), 'color',colorarray(i,:))
    % hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'o', 'color',colorarray(i,:))
    hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'bo','linewidth',2,'markersize',7);
    
end
xlabel(['Average Branch Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'FontSize',13);
ylabel('Cell Speed','FontSize',13);
title(['Branchness vs Speed, ',...
    ' Correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean', ...
    Group_Pool_Travel_Speed'),'%1.2f')],'FontSize',13);
% axis([0 10 0 24]);
set(gca,'fontsize',13);
saveas(h7,[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.fig']);
saveas(h7,[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.tif']);
print(h7,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.eps']);


h8 = figure(8);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16)
    plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
end
xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'Fontsize',13);
ylabel('Cell Vim mean Int','Fontsize',13);
title(['Branch number vs Vim, correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
set(gca,'fontsize',13);
saveas(h8,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.fig']);
saveas(h8,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.tif']);
print(h8,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.eps']);


h8_axis = axis;

h18 = figure(18);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    %     text(Group_Pool_thresholded_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
end
corr_b_v_v = corrcoef(Group_Pool_whole_cell_vif_mean_intensity',Group_Pool_thresholded_branch_number_mean');
display(['Corr for Branchness vs Vim: ',num2str(corr_b_v_v(1,2))]);

xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_thresholded_branch_number_mean))],'Fontsize',13);
ylabel('Cell Vim mean Int','Fontsize',13);
title({['Branch number vs Vim, for branches larger than ',num2str(T_branchsize),' pixels'],['duration longer than ',num2str(T_branchduration),' frames'],...
    ['Correlation:',num2str(corr_b_v_v(1,2))]},'Fontsize',13);
axis(h8_axis);
set(gca,'fontsize',13);
saveas(h18,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.fig']);
saveas(h18,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.tif']);
print(h18,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.eps']);



%%
try
    %%
    h26 = figure(26);hold off;
    % plot(Group_Pool_fila_branch_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
        % text(Group_Pool_fila_branch_orientation_pool_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
        hold on;
        % plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
        
    end
    % axis([0.55, 1.0, 350 950]);
    xlabel('Filament Orientation along Branch, Standard Deviation','Fontsize',13);
    ylabel('Cell Vim mean Int','Fontsize',13);
    title({'Filament Orientation Scatterness vs Vim Mean Int',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
        Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f'), ...
        ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
    set(gca,'fontsize',13);
    saveas(h26,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.fig']);
    saveas(h26,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.tif']);
    
    print(h26,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.eps']);
    
    
    %%
    
    h36 = figure(36);hold off;
    % plot(Group_Pool_fila_branch_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_branch_number_mean)
        % text(Group_Pool_fila_branch_orientation_pool_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
        hold on;
        %         plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_branch_number_mean(i), 'bo','linewidth',2,'markersize',7);
    end
    xlabel('Filament Orientation along Branch, Standard Deviation','Fontsize',13);
    ylabel('Branch Mean Number ','Fontsize',13);
    title({'Filament Orient along Branch Scatterness vs Branch Mean Num',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))', ...
        Group_Pool_branch_number_mean(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f'),...
        ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))))]},'Fontsize',13);
    axis([0.6 0.9 0 10]);
    set(gca,'fontsize',13);

    saveas(h36,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.fig']);
    saveas(h36,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.tif']);
    print(h36,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.eps']);
    
    
end

try
    h46 = figure(46);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    
    title('Filament Orientation along Branch Orientaion Distribution','Fontsize',13,'fontname','Helvetica');
    set(gca,'fontsize',13);
    saveas(h46,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.fig']);
    saveas(h46,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.tif']);
    
    print(h46,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.eps']);
    
    
    
    h56 = figure(56);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool_fast,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    title('Filament Orientation along Branch Orientaion Distribution(fast cell)','Fontsize',13,'fontname','Helvetica');
    set(gca,'fontsize',13);
    saveas(h56,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.fig']);
    saveas(h56,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.tif']);
    print(h56,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.eps']);
    
    
    h66 = figure(66);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool_slow,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    set(gca,'fontsize',13);
    title('Filament Orientation along Branch Orientaion Distribution(slow cell)','Fontsize',13,'fontname','Helvetica');
    
    saveas(h66,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.fig']);
    saveas(h66,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.tif']);
    print(h66,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.eps']);
    
    
    %%
    
    %%
    h27 = figure(27);hold off;
    % plot(Group_Pool_fila_trajectory_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
        % text(Group_Pool_fila_trajectory_orientation_pool_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
        hold on;
        % plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
        
    end
    % axis([0.55, 1.0, 350 950]);
    xlabel('Filament Orientation along Path, Standard Deviation','Fontsize',13);
    ylabel('Cell Vim mean Int','Fontsize',13);
    title({'Filament Orientation along Path Scatterness vs Vim Mean Int',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
        Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f'),...
        ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
    set(gca,'fontsize',13);
    saveas(h27,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.fig']);
    saveas(h27,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.tif']);
    print(h27,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.eps']);
    
    
    %%
end

try
    h37 = figure(37);hold off;
    % plot(Group_Pool_fila_trajectory_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_branch_number_mean)
        % text(Group_Pool_fila_trajectory_orientation_pool_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
        hold on;
        %         plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_branch_number_mean(i),...
            'bo','linewidth',2,'markersize',7);
    end
    xlabel('Filament Orientation along Cell Movement, Standard Deviation','Fontsize',13);
    ylabel('Branch Mean Number ','Fontsize',13);
    title({'Filament Orient along Cell Path Scatterness vs Branch Mean Num',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))', ...
        Group_Pool_branch_number_mean(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f') ,...
        ', Sample Size:',num2str(numel(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))))]},'Fontsize',13);
    set(gca,'fontsize',13);
    saveas(h37,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.fig']);
    saveas(h37,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.tif']);
    
    print(h37,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.eps']);
    
    
end

try
    h57 = figure(57);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    
    title('Filament Orientation along Cell Movement Distribution','Fontsize',13,'fontname','Helvetica');
    
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.fig']);
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.tif']);
    
    print(h57,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.eps']);
    
    h57 = figure(57);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool_fast,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    title('Filament Orientation along Cell Movement Distribution(fast cell)','Fontsize',13,'fontname','Helvetica');
    
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.fig']);
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.tif']);
    
    print(h57,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.eps']);
    
    
    h67 = figure(67);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool_slow,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    title('Filament Orientation along Cell Movement Distribution(slow cells)','Fontsize',13,'fontname','Helvetica');
    set(gca,'fontsize',13);
    saveas(h67,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.fig']);
    saveas(h67,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.tif']);
    print(h67,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.eps']);
    
end
%%
h28 = figure(28);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    % plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
    
end
% axis([0.55, 1.0, 350 950]);
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
ylabel('Cell Vim mean Int','Fontsize',13);
title({'Branch Orientation Scatterness vs Vim',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
    Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f') ,...
    ', Sample Size:',num2str(numel(Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_whole_cell_vif_mean_intensity)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h28,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.fig']);
saveas(h28,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.tif']);
print(h28,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.eps']);


%%

h38 = figure(38);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_branch_number_mean)
    % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
ylabel('Branch Mean Number ','Fontsize',13);
title({'Branch Orientation Scatterness vs Branch Mean Number',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_branch_number_mean))', ...
    Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f') ,...
    ', Sample Size:',num2str(numel(Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_branch_cellmovement_std))))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h38,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.fig']);
saveas(h38,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.tif']);
print(h38,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.eps']);

h48 = figure(48);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : numel(Group_Pool_whole_cell_vif_mean_intensity)
    text(Group_Pool_Travel_Speed(i)-0.02, ...
        Group_Pool_whole_cell_vif_mean_intensity(i)+0.02, ...
        Identifier_cell{i}, 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Cell Speed','Fontsize',13);
ylabel('Vim Level ','Fontsize',13);
title({'Cell Speed vs Vim Level',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_Travel_Speed', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f'),...
    ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h48,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.fig']);
saveas(h48,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.tif']);
print(h48,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.eps']);



h58 = figure(58);hold off;

% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_Travel_Speed(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+0.2, num2str(i), 'color',colorarray(i,:));
    if(Group_Pool_Cell_Marked_Frame_Number(i)>15 && Group_Pool_CompletedFrame_last(i)>60)
    hold on;
    %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
    end
end
xlabel('Cell Speed','Fontsize',13);
ylabel('Vim Level ','Fontsize',13);
title({'Cell Speed vs Vim Level (cell with more frames, with late hours)',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)', ...
    Group_Pool_whole_cell_vif_mean_intensity(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)'),'%1.2f'),...
    ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h58,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.fig']);
saveas(h58,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.tif']);
print(h58,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.eps']);



%%
%
h9 = figure(9);
plot(Group_Pool_branch_number_mean_pat, Group_Pool_branch_vif_mean_intensity,'.');
xlabel('Branch Mean Number','fontsize',13);
ylabel('Branch Vim mean Int','fontsize',13);
title('Branch number vs Vim','fontsize',13);

set(gca,'fontsize',13);
saveas(h9,[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.fig']);
saveas(h9,[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.tif']);
print(h9,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.eps']);


h10 = figure(10);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vim_totalamount_mean)
    %     text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
    plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), 'bo','linewidth',2,'markersize',7);
    
end
xlabel('Branch Mean Number','fontsize',13);
ylabel('Cell Vim Total Int','fontsize',13);
title({['Branch number vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')'],...
    ['Correlation: ',num2str(corr(Group_Pool_branch_number_mean',Group_Pool_whole_cell_vim_totalamount_mean'))]},'fontsize',13);
set(gca,'fontsize',13);

saveas(h10,[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.fig']);
saveas(h10,[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.tif']);
print(h10,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.eps']);


% h20 = figure(20);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% hold on;
% plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% end
% xlabel('Branch Mean Number','Fontsize',13);
% ylabel('Cell Vim Mean Int','Fontsize',13);
% title(['Branch number vs Vim Mean,  correlation: ',...
%     num2str(corr(Group_Pool_branch_number_mean', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% saveas(h20,[Group_ROOT_DIR,filesep,'Branchness_vs_VimMean.fig']);
% saveas(h20,[Group_ROOT_DIR,filesep,'Branchness_vs_VimMean.tif']);
%
% h30 = figure(30);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% hold on;
% plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% end
% xlabel('Cell Speed','Fontsize',13);
% ylabel('Cell Vim Mean Int','Fontsize',13);
% title(['Cell Speed vs Vim Mean,  correlation: ',...
%     num2str(corr(Group_Pool_Travel_Speed', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% saveas(h30,[Group_ROOT_DIR,filesep,'Speed_vs_VimMean.fig']);
% saveas(h30,[Group_ROOT_DIR,filesep,'Speed_vs_VimMean.tif']);


h11 = figure(11);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_tracked_branch_d_frames)
    %     text(Group_Pool_tracked_branch_d_frames(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
    plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Tracked Branchs per Frame','fontsize',13);
ylabel('Cell Vim Total Int','fontsize',13);
title({['Tracked Branches per Frame vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')'],...
    ['Correlation: ',num2str(corr(Group_Pool_tracked_branch_d_frames',Group_Pool_whole_cell_vim_totalamount_mean'))]},'fontsize',13);
set(gca,'fontsize',13);

saveas(h11,[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.fig']);
saveas(h11,[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.tif']);
print(h11,'-depsc',[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.eps']);


h12 = figure(12);
plot(Group_Pool_whole_cell_vif_mean_intensity_pat, Group_Pool_branch_duration_array, '.');
ylabel('Branch Duration','fontsize',13);
xlabel('Cell Vim mean Int','fontsize',13);
title('Cell Vim vs Branch Duration','fontsize',13);
set(gca,'fontsize',13);

saveas(h12,[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.fig']);
saveas(h12,[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.tif']);
print(h12,'-depsc',[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.eps']);


h13 = figure(13);
plot(Group_Pool_branch_duration_array(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)), ...
    Group_Pool_branch_vif_mean_intensity(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)),'.');
xlabel('Branch Duration','fontsize',13);
ylabel('Branch Vim mean Int','fontsize',13);
title({['Branch: Duration vs Vim, for branches larger than ',num2str(T_branchsize),' pixels,'],['duration longer than ',num2str(T_branchduration),' frames']},'fontsize',13);
set(gca,'fontsize',13);
saveas(h13,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.fig']);
saveas(h13,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.tif']);
print(h13,'-depsc',[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.eps']);



%%
branch_analysis_plotting_VimNms;
branch_analysis_plotting_VimTotalNms;
branch_analysis_plotting_VimFilaDen;
branch_analysis_plotting_VimFilaTotal;

save([Group_ROOT_DIR,filesep,'branch_analysis_group_results.mat']);

