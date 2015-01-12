function BA_output_ML_cell = branch_analysis_moiveList_results_gather(ML, T_branchsize, T_branchduration, Group_ROOT_DIR)
% function to gather information of branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML:     The movieList object to be analyzed

%% take the input default
if(nargin<2)
    T_branchsize=200;
end
if(nargin<3)
    T_branchduration=2;
end
if(nargin<4)
    Group_ROOT_DIR=[];
end

%% correct 0 input for these parameters
if(T_branchsize<=0)
    T_branchsize=200;
end

if(T_branchduration<=0)
    T_branchduration=2;
end

%% settings
ML_ROOT_DIR = ML.outputDirectory_;

if(isempty(Group_ROOT_DIR))
    Group_ROOT_DIR = ML_ROOT_DIR;
end

% the number of movies
movieNumber =  length(ML.movieDataFile_);

% load one movie by one movie
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
            outputPath = [ROOT_DIR,'\BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];
            
            % check if this folder exist
            if(exist(outputPath,'dir'))
                % if it exist, try to do the branch analysis
                try
                    if(exist([outputPath,'\branch_analysis_results.mat'],'file'))
                        load([outputPath,'\branch_analysis_results.mat'],'BA_output');
                    else
                        continue;
                    end
                    
                    display(['Found: iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                    
                    % save
                    BA_output_ML_cell{1, iM}{iChannel, iCell} = BA_output;
                end
            end
        end
    end
end

%% Gather them into pools

% initialize the pools

Group_Pool_Travel_Length = [];
Group_Pool_Travel_Distance = [];
Group_Pool_Travel_Speed = [];
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

iAllCell = 0;


for iM  = 1 : movieNumber
    
    for iChannel = 1 : nChannel
        
        for iCell = 1 :  10
            
            try
                BA_output = BA_output_ML_cell{1, iM}{iChannel, iCell};
            catch
                %  there is no such cell in this channel
                continue;
            end
            
            % if the result for this channel this cell is not empty
            % pool them up
            if(~isempty(BA_output))
                Group_Pool_branch_size_mean = [Group_Pool_branch_size_mean BA_output.branch_mean_size];
                Group_Pool_whole_cell_size_mean= [Group_Pool_whole_cell_size_mean BA_output.whole_cell_size_mean];
                Group_Pool_branch_number_mean_pat = [Group_Pool_branch_number_mean_pat repmat(BA_output.branch_number_mean, [1, length(BA_output.branch_vif_mean_intensity)])];
                Group_Pool_whole_cell_vif_mean_intensity_pat = [Group_Pool_whole_cell_vif_mean_intensity_pat repmat(BA_output.whole_cell_vif_mean_intensity, [1, length(BA_output.branch_vif_mean_intensity)])];
                Group_Pool_whole_cell_vim_totalamount_mean = [Group_Pool_whole_cell_vim_totalamount_mean BA_output.whole_cell_vim_totalamount_mean];
                
                
                Thresholded_branch_number_mean = sum(BA_output.branch_duration_array(find(BA_output.branch_mean_size>T_branchsize & BA_output.branch_duration_array>T_branchduration)))...
                    ./ BA_output.cell_marked_frame_number;
                Group_Pool_thresholded_branch_number_mean = [Group_Pool_thresholded_branch_number_mean Thresholded_branch_number_mean];
                
                Group_Pool_Travel_Length = [Group_Pool_Travel_Length BA_output.cell_travel_length];
                Group_Pool_Travel_Distance = [Group_Pool_Travel_Distance BA_output.cell_travel_distance];
                Group_Pool_Travel_Speed = [Group_Pool_Travel_Speed BA_output.cell_travel_speed];
                Group_Pool_Cell_Marked_Frame_Number = [Group_Pool_Cell_Marked_Frame_Number BA_output.cell_marked_frame_number];
                Group_Pool_branch_number_tracked = [Group_Pool_branch_number_tracked BA_output.branch_number_tracked];
                Group_Pool_branch_number_max = [Group_Pool_branch_number_max BA_output.branch_number_max];
                Group_Pool_branch_number_mean = [Group_Pool_branch_number_mean BA_output.branch_number_mean];
                Group_Pool_branch_duration_array = [Group_Pool_branch_duration_array  BA_output.branch_duration_array];
                Group_Pool_branch_duration_mean= [Group_Pool_branch_duration_mean  BA_output.branch_duration_mean];
                Group_Pool_branch_vif_mean_intensity= [Group_Pool_branch_vif_mean_intensity BA_output.branch_vif_mean_intensity];
                Group_Pool_protrusion_vif_mean_intensity= [Group_Pool_protrusion_vif_mean_intensity BA_output.protrusion_vif_mean_intensity];
                Group_Pool_retraction_vif_mean_intensity= [Group_Pool_retraction_vif_mean_intensity BA_output.retraction_vif_mean_intensity];
                Group_Pool_whole_cell_vif_mean_intensity= [Group_Pool_whole_cell_vif_mean_intensity BA_output.whole_cell_vif_mean_intensity];
                Group_Pool_branch_cellmovement_std= [Group_Pool_branch_cellmovement_std BA_output.branch_cellmovement_std];
                
                iAllCell = iAllCell + 1;
                
            end
        end        
    end
end

Group_Pool_tracked_branch_d_frames = Group_Pool_branch_number_tracked./Group_Pool_Cell_Marked_Frame_Number;

save([Group_ROOT_DIR,'\movieList_BA_results_gathered.mat']);

