function [BA_output, f_display] = branch_analysis_marked_cell_GMST(MD, iChannel, iCell, half_size,min_branch_size_Threshold,filament_stat_flag,figure_flag,min_filament_length)
% function to do branch analysis for marked cells
% Liya Ding, Feb, 2014
%
% Input:
%   MD:     The movieData object loaded before running this function
%   iChannel: which channel is the membrane, so the other is the VIM by default
%   iCell: the index of marked cell (for single cell segmentation)
%
% Output:
%   BA_output: is a struct with these fields:
%           cell_travel_length:   the cell center travelled along the curved path, this is the length of the path
%           cell_travel_distance: the begin to end position distance
%           cell_travel_speed: the speed of migration, here as along the path
%           cell_marked_frame_number: the number of frames used.
%           branch_number_tracked: the number of branches tracked
%           branch_number_max: each frame has a number of branches, take the maximum among all the frames
%           branch_number_mean: the average branches number among all the frames
%           branch_duration_array: the branch by branch duration value
%           branch_duration_mean: the average duration value
%           branch_vif_mean_intensity: the branch by branch, (average within the same branch) vif intensity value
%           protrusion_vif_mean_intensity: the average vim intensity in the protrusion(the red) region, among all frames.
%           retraction_vif_mean_intensity: the average vim intensity in the retraction(the green) region, among all frames.

f_display=[];

if(nargin<4)
    half_size=150;
end

if(nargin<5)
    min_branch_size_Threshold=100;
end

if(nargin<6)
    filament_stat_flag=0;
end


if(nargin<7)
    figure_flag=0;
end

if(nargin<8)
    min_filament_length=0;
end

display(['iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);

display_msg_flag=0;
package_process_ind_script;

%initialize an empty output
BA_output=[];


if(iChannel==1)
    VIF_channel = 2;
else
    VIF_channel = 1;
end

if numel(MD.channels_) < VIF_channel
    VIF_channel = numel(MD.channels_);
end

ROOT_DIR = MD.outputDirectory_;

FilamentAnalysisPackage_complete_frames_file_name = [ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),filesep,'completedFrames.mat'];
SegmentationPackage_complete_frames_file_name = [ROOT_DIR,filesep,'SegmentationPackage',filesep,'completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),filesep,'completedFrames.mat'];

if(exist(FilamentAnalysisPackage_complete_frames_file_name,'file'))
    PackageName = 'FilamentAnalysisPackage';
end

if(exist(SegmentationPackage_complete_frames_file_name,'file'))
    PackageName = 'SegmentationPackage';
end


load([ROOT_DIR,filesep,PackageName,filesep,'completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),filesep,'completedFrames.mat'],'isCompleted');
truthPath = [ROOT_DIR,filesep,PackageName,filesep,'OnlyFixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
old_truthPath = [ROOT_DIR,filesep,PackageName,filesep,'FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
outputPath = [ROOT_DIR,filesep,'BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];

if(sum(isCompleted)==0)
    return;
end

fileExistance = 0*isCompleted;

for iFrame = 1 : MD.nFrames_
    if(exist([old_truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file') ...
            || exist([truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file'))
        fileExistance(iFrame)=1;
    end
end

isCompleted = (isCompleted.*fileExistance)>0;

% only use the longest continuously marked sequence
isCompleted_longest = keep_largest_area(isCompleted);

if(sum(isCompleted -isCompleted_longest)>15)
    display('Please check, there is more sequence available');
end

isCompleted = isCompleted_longest;

if(sum(isCompleted)<3)
    return;
end

if(~exist(outputPath,'dir'))
    mkdir(outputPath);
end

CompletedFrame = find(isCompleted>0);
CompletedFrame = CompletedFrame(:);
CompletedFrame = CompletedFrame';


% keep the frame numbers
BA_output.CompletedFrame = CompletedFrame;

% these index is sure to be continuous
FirstFrame = min(CompletedFrame);
nCompleteFrame = length(CompletedFrame);

raw_mask_cell = cell(1,nCompleteFrame);
smoothed_mask_cell = cell(1,nCompleteFrame);
current_seg_cell_st = cell(1,nCompleteFrame);
orienation_map_filtered_cell_st = cell(1,nCompleteFrame);
current_seg_cell_gm = cell(1,nCompleteFrame);
orienation_map_filtered_cell_gm = cell(1,nCompleteFrame);
region_orientation_cell = cell(1,nCompleteFrame);
current_model_cell = cell(1,nCompleteFrame);


min_y = Inf;
min_x = Inf;
max_y = 1;
max_x = 1;
Channel_FilesNames = MD.channels_(VIF_channel).getImageFileNames(1:MD.nFrames_);

filename_short_strs = uncommon_str_takeout(Channel_FilesNames);


new_CompletedFrame = [];

for iFrame = CompletedFrame
    display(['iChannel:', num2str(iChannel),', iCell:', num2str(iCell),', iFrame:',num2str(iFrame) ]);
    iCompleteFrame = iFrame - FirstFrame +1;
    
    %% filament segmentation part
    current_seg_gm=[];
    orienation_map_filtered_gm=[];
    % default as empty in filament segmentation
    current_seg_st=[];
    orienation_map_filtered_st=[];
    nms=[];
    
    % GM
    % if filament stat is available and requested
    if(exist(GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation_GM',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'filament_seg_',filename_short_strs{iFrame},'.mat']),'file') && filament_stat_flag>0)
        load( GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation_GM',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'filament_seg_',...
            filename_short_strs{iFrame},'.mat']),'current_seg_orientation','current_model');
        
        current_seg_gm = ~isnan(current_seg_orientation);
        orienation_map_filtered_gm = current_seg_orientation;
        [length_limited_model,length_limited_seg] = ...
            filament_model_length_limit(current_model,current_seg_gm,min_filament_length);
        current_seg_gm =  ~isnan(length_limited_seg);
        
    else
        if(exist(GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation_GM',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'steerable_vote_',filename_short_strs{iFrame},'.mat']),'file') && filament_stat_flag>0)
            load( GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation_GM',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'steerable_vote_',...
                filename_short_strs{iFrame},'.mat']),'current_seg_orientation','current_model');
            current_seg_gm = ~isnan(current_seg_orientation);
            orienation_map_filtered_gm = current_seg_orientation;
          [length_limited_model,length_limited_seg] = ...
            filament_model_length_limit(current_model,current_seg_gm,min_filament_length);
        current_seg_gm =  ~isnan(length_limited_seg);
          
        end
    end
    
    % ST
    % if filament stat is available and requested
    if(exist(GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'filament_seg_',filename_short_strs{iFrame},'.mat']),'file') && filament_stat_flag>0)
        load( GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'filament_seg_',...
            filename_short_strs{iFrame},'.mat']),'current_seg_orientation','current_model');
        current_seg_st = ~isnan(current_seg_orientation);
        orienation_map_filtered_st = current_seg_orientation;        
        
    else
        if(exist(GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'steerable_vote_',filename_short_strs{iFrame},'.mat']),'file') && filament_stat_flag>0)
            load( GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'steerable_vote_',...
                filename_short_strs{iFrame},'.mat']),'current_seg_orientation','current_model');
            current_seg_st = ~isnan(current_seg_orientation);
            orienation_map_filtered_st = current_seg_orientation;
            
        end
    end
    
    
   
    % if filament stat is available and requested
    if(exist(GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'SteerableFiltering',filesep,'Channel',num2str(VIF_channel),filesep,'steerable_',filename_short_strs{iFrame},'.mat']),'file') && filament_stat_flag>0)
        load( GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'SteerableFiltering',filesep,'Channel',num2str(VIF_channel),filesep,'steerable_',...
            filename_short_strs{iFrame},'.mat']),'nms', 'MAX_st_res','scaleMap');
    end
    
         
    %% Load the mask and smooth the cell mask
    if(exist([truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file'))
        current_mask = imread([truthPath,filesep,'mask_',num2str(iFrame),'.tif']);
    else
        if(exist([old_truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file'))
            current_mask = imread([old_truthPath,filesep,'mask_',num2str(iFrame),'.tif']);
        else
            disp('no mask found');
            break;
        end
    end
    
    current_seg_cell_st{1,iCompleteFrame} = current_seg_st.*double(current_mask);    
    orienation_map_filtered_cell_st{1,iCompleteFrame} = orienation_map_filtered_st.*double(current_mask);
    current_seg_cell_gm{1,iCompleteFrame} = current_seg_gm.*double(current_mask); 
    orienation_map_filtered_cell_gm{1,iCompleteFrame} = orienation_map_filtered_gm.*double(current_mask);
    
    
    %%%%% get the st amount and gm amount and the ratio
 %   if(mod(iFrame,4)==1)        
        st_sum(iCompleteFrame) = sum(sum(current_seg_cell_st{1,iCompleteFrame}));
        gm_sum(iCompleteFrame)  = sum(sum(current_seg_cell_gm{1,iCompleteFrame}));
        gm_in_st_sum(iCompleteFrame) = sum(sum(current_seg_cell_st{1,iCompleteFrame}.*current_seg_cell_gm{1,iCompleteFrame}));
        gm_st_ratio(iCompleteFrame) = gm_sum/st_sum;
        gminst_st_ratio(iCompleteFrame) = gm_in_st_sum/st_sum;        
  %  end
            
    current_mask([ 1 end], :)=0;
    current_mask(:, [ 1 end])=0;
    current_mask = keep_largest_area(current_mask);
    raw_mask_cell{1,iCompleteFrame} = current_mask>0;
    [indy,indx]= find(current_mask>0);
    min_y = min(min_y,min(indy));
    min_x = min(min_x,min(indx));
    max_y = max(max_y,max(indy));
    max_x = max(max_x,max(indx));
    
    B = bwboundaries(current_mask);
    
    Y = B{1}(:,1);
    X = B{1}(:,2);
    
    H = fspecial('gaussian',2*half_size+1,half_size/4);
    H = double(H(:,half_size+1));
    H = H./(sum(H));
    
    % smooth in a circular fashion, since this is an enclosed curve
    X = imfilter(X, H,'circular','same');
    Y = imfilter(Y, H,'circular','same');
    X = imfilter(X, H,'circular','same');
    Y = imfilter(Y, H,'circular','same');
    
    % then rebuild the mask
    smoothed_current_mask = roipoly(current_mask,X,Y);
    
    % if there is nothing in cell mask, break here
    if(sum(sum(smoothed_current_mask))==0)
        disp('empty mask');
        break;
    end
    
    % keep the smoothed mask in the cell-structure
    smoothed_current_mask = keep_largest_area(smoothed_current_mask);
    
    % make it logic to save memory
    smoothed_mask_cell{1,iCompleteFrame} = (smoothed_current_mask>0);
    
    % keep track if stopped before the loop ends
    new_CompletedFrame = [new_CompletedFrame iFrame];
end

% correcting the early stops
CompletedFrame = new_CompletedFrame;
nCompleteFrame = length(CompletedFrame);


% for iCompleteFrame = 1 : nCompleteFrame
%     imwrite(smoothed_mask_cell{1,iCompleteFrame}, [outputPath,filesep,'smoothed_marked_mask_',num2str(iCompleteFrame),'.tif']);
% end

color_array = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 1 0; 0 1 1; rand(494,3)];
region_branch_label_cell = cell(1,500);
label_skel_cell = cell(1,500);
branch_leaf_flag_cell= cell(1,500);
% region_branch_wo_inner_label_cell= cell(1,500);
% label_skel__wo_inner_cell{iCompleteFrame}= cell(1,500);

for iCompleteFrame = 1 : nCompleteFrame
    display(['iChannel:', num2str(iChannel),', iCell:', num2str(iCell),', iCompleteFrame:',num2str(iCompleteFrame) ]);
      
    smoothed_current_mask = smoothed_mask_cell{1,iCompleteFrame};
    
    try
       % change the way of determine the cell center as the Hunter's idea
        % in his paper, but briefly in 2D
             
        % invert
        inverted_Distmap = (bwdist(1-smoothed_current_mask));
        % get the maximum points
        max_invertdist = inverted_Distmap==max(max(inverted_Distmap));
        
        % in case there are more than one of these max centers
        max_invertdist = keep_largest_area(max_invertdist);
        
        % in each maximum, there could be more than one point, to make the
        % code easier for its own flow
        C_xy = regionprops(max_invertdist,'Centroid');
        center_x(iCompleteFrame) = C_xy.Centroid(1);
        center_y(iCompleteFrame) = C_xy.Centroid(2);
%         
    catch
        disp('error in cell mask centroid');
        break;
    end
    
end

BA_output.cell_travel_length = ...
    sum(sqrt((center_x(2:end)-center_x(1:end-1)).^2 + (center_y(2:end)-center_y(1:end-1)).^2 ));

BA_output.cell_travel_distance = ...
    (sqrt((center_x(end)-center_x(1)).^2 + (center_y(end)-center_y(1)).^2 ));

BA_output.cell_travel_speed = ...
    BA_output.cell_travel_length./nCompleteFrame;

BA_output.cell_marked_frame_number = ...
    nCompleteFrame;


      BA_output.st_sum = st_sum;
      BA_output.gm_sum = gm_sum;
       BA_output.gm_in_st_sum = gm_in_st_sum;
      BA_output.gm_st_ratio =gm_st_ratio ;
       BA_output.gminst_st_ratio = gminst_st_ratio;
     
   

%% get cell centroid trajectory direction
center_x(iCompleteFrame) = C_xy.Centroid(1);
center_y(iCompleteFrame) = C_xy.Centroid(2);

% take 1/5 of the smoothing size of the cell boundary
trajectory_smooth_size = round(half_size/5);
H = fspecial('gaussian',2*trajectory_smooth_size+1,trajectory_smooth_size/4);
H = double(H(:,trajectory_smooth_size+1));
H = H./(sum(H));

smooth_center_x =  imfilter(center_x,H','replicate','same');
smooth_center_y =  imfilter(center_y,H','replicate','same');

%%
BA_output.center_x = center_x;
BA_output.center_y = center_y;
BA_output.smooth_center_x = smooth_center_x;
BA_output.smooth_center_y = smooth_center_y;
BA_output.speed_marked_frames = sqrt((center_x(2:end)-center_x(1:end-1)).^2 + (center_y(2:end)-center_y(1:end-1)).^2 );
BA_output.smoothed_speed_marked_frames = sqrt((smooth_center_x(2:end)-smooth_center_x(1:end-1)).^2 + (smooth_center_y(2:end)-smooth_center_y(1:end-1)).^2 );

% persistency calculation



if(length(smooth_center_x)>3)
    
    for iF = 1 : length(smooth_center_x)-2
        persistency3(iF) = sqrt((smooth_center_x(iF+2)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+2)-smooth_center_y(iF)).^2 )...
            /sum(sqrt((smooth_center_x(iF+1:iF+2)-smooth_center_x(iF:iF+1)).^2 + (smooth_center_y(iF+1:iF+2)-smooth_center_y(iF:iF+1)).^2));
    end
    BA_output.persistency3=persistency3;
    
end

if(length(smooth_center_x)>4)    
    for iF = 1 : length(smooth_center_x)-3
        persistency4(iF) = sqrt((smooth_center_x(iF+3)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+3)-smooth_center_y(iF)).^2 )...
            /sum(sqrt((smooth_center_x(iF+1:iF+3)-smooth_center_x(iF:iF+2)).^2 + (smooth_center_y(iF+1:iF+3)-smooth_center_y(iF:iF+2)).^2));
    end
    BA_output.persistency4=persistency4;    
end

if(length(smooth_center_x)>5)
    
    for iF = 1 : length(smooth_center_x)-4
        persistency5(iF) = sqrt((smooth_center_x(iF+4)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+4)-smooth_center_y(iF)).^2 )...
            /sum(sqrt((smooth_center_x(iF+1:iF+4)-smooth_center_x(iF:iF+3)).^2 + (smooth_center_y(iF+1:iF+4)-smooth_center_y(iF:iF+3)).^2));
    end
    BA_output.persistency5=persistency5;
    
end

if(length(smooth_center_x)>6)
    
    for iF = 1 : length(smooth_center_x)-5
        persistency6(iF) = sqrt((smooth_center_x(iF+5)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+5)-smooth_center_y(iF)).^2 )...
            /sum(sqrt((smooth_center_x(iF+1:iF+5)-smooth_center_x(iF:iF+4)).^2 + (smooth_center_y(iF+1:iF+5)-smooth_center_y(iF:iF+4)).^2));
    end
    BA_output.persistency6=persistency6;
    
end

% pace of 15 since the min length is 20
for iF = 1 : length(smooth_center_x)-14
   persistency15(iF) = sqrt((smooth_center_x(iF+14)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+14)-smooth_center_y(iF)).^2 )...
       /sum(sqrt((smooth_center_x(iF+1:iF+14)-smooth_center_x(iF:iF+13)).^2 + (smooth_center_y(iF+1:iF+14)-smooth_center_y(iF:iF+13)).^2));
end
BA_output.persistency15=persistency15;


if(length(smooth_center_x)>20)
    for iF = 1 : length(smooth_center_x)-19
        persistency20(iF) = sqrt((smooth_center_x(iF+19)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+19)-smooth_center_y(iF)).^2 )...
            /sum(sqrt((smooth_center_x(iF+1:iF+19)-smooth_center_x(iF:iF+18)).^2 + (smooth_center_y(iF+1:iF+19)-smooth_center_y(iF:iF+18)).^2));
    end
    BA_output.persistency20=persistency20;
end



if(length(smooth_center_x)>30)
    
    for iF = 1 : length(smooth_center_x)-29
        persistency30(iF) = sqrt((smooth_center_x(iF+29)-smooth_center_x(iF)).^2 + (smooth_center_y(iF+29)-smooth_center_y(iF)).^2 )...
            /sum(sqrt((smooth_center_x(iF+1:iF+29)-smooth_center_x(iF:iF+28)).^2 + (smooth_center_y(iF+1:iF+29)-smooth_center_y(iF:iF+28)).^2));
    end
    BA_output.persistency30=persistency30;
    
end



%%

save([outputPath,filesep,'branch_analysis_results_gmst.mat'],'BA_output');

