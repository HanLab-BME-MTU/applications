function [cell_size_mean, total_vim_mean] = cell_size_after_branch_analysis_marked_cell(MD, iChannel, iCell)
% function to add the cell size after doing branch analysis for marked cells
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


package_process_ind_script;

if(iChannel==1)
    VIF_channel = 2;
else
    VIF_channel = 1;
end

ROOT_DIR = MD.outputDirectory_;

if(exist([ROOT_DIR,'\FilamentAnalysisPackage\refined_masks\'],'dir'))
    PackageName = 'FilamentAnalysisPackage';
else
    PackageName = 'SegmentationPackage';
end

load([ROOT_DIR,'\',PackageName,'\completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),'\completedFrames.mat'],'isCompleted');
truthPath = [ROOT_DIR,'\',PackageName,'\FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
outputPath = [ROOT_DIR,'\BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];

if(~exist(outputPath,'dir'))
    mkdir(outputPath);
end


% only use the longest continuously marked sequence
isCompleted = keep_largest_area(isCompleted);

CompletedFrame = find(isCompleted>0);
CompletedFrame = CompletedFrame(:);
CompletedFrame = CompletedFrame';

% these index is sure to be continuous
FirstFrame = min(CompletedFrame);
nCompleteFrame = length(CompletedFrame);

raw_mask_cell = cell(1,nCompleteFrame);
smoothed_mask_cell = cell(1,nCompleteFrame);
current_seg_cell= cell(1,nCompleteFrame);
orienation_map_filtered_cell= cell(1,nCompleteFrame);

min_y = Inf;
min_x = Inf;
max_y = 1;
max_x = 1;

for iFrame = CompletedFrame
    iCompleteFrame = iFrame - FirstFrame +1;
    
    current_mask = imread([truthPath,'/mask_',num2str(iFrame),'.tif']);
    current_mask([ 1 end], :)=0;
    current_mask(:, [ 1 end])=0;
    current_mask = keep_largest_area(current_mask);
    raw_mask_cell{1,iCompleteFrame} = current_mask>0;
    [indy,indx]= find(current_mask>0);
    min_y = min(min_y,min(indy));
    min_x = min(min_x,min(indx));
    max_y = max(max_y,max(indy));
    max_x = max(max_x,max(indx));
    
%     load([ROOT_DIR,'/FilamentAnalysisPackage/FilamentSegmentation/Channel2/DataOutput/steerable_vote_',...
%         num2str(iFrame),'.mat'],'current_seg','orienation_map_filtered');
%     
     current_seg_cell{1,iCompleteFrame} = [];
%     orienation_map_filtered_cell{1,iCompleteFrame} = orienation_map_filtered;
%     
%     %     imwrite(current_mask, ['nonboarder_mask_',num2str(iFrame),'.tif']);
%     
    B = bwboundaries(current_mask);
    
    Y = B{1}(:,1);
    X = B{1}(:,2);
    
    H = fspecial('gaussian',121,15);
    H = double(H(:,61));
    H = H./(sum(H));
    
    X = imfilter(X, H,'replicate','same');
    Y = imfilter(Y, H,'replicate','same');
    X = imfilter(X, H,'replicate','same');
    Y = imfilter(Y, H,'replicate','same');
    
    smoothed_current_mask = roipoly(current_mask,X,Y);
    
    % make it logic to save memory
    smoothed_mask_cell{1,iCompleteFrame} = (smoothed_current_mask>0);
    
end


cell_size_pool = [];
cell_vimtotal_pool = [];

for iCompleteFrame = 1 : nCompleteFrame
    %     current_seg = current_seg_cell{1,iFrame};
    iFrame = iCompleteFrame+FirstFrame-1;
    
    smoothed_current_mask = smoothed_mask_cell{1,iCompleteFrame};
    
    current_VIF_image = MD.channels_(VIF_channel). loadImage(iFrame);
    
    cell_size_pool = [cell_size_pool; sum(sum(smoothed_current_mask>0))];
    cell_vimtotal_pool = [cell_vimtotal_pool; sum(current_VIF_image(smoothed_current_mask>0))];
        
    end

load([outputPath,'\branch_analysis_results.mat'],'BA_output');

BA_output.whole_cell_vim_totalamount_mean = mean(cell_vimtotal_pool);

BA_output.whole_cell_size_mean = mean(cell_size_pool);

save([outputPath,'\branch_analysis_results.mat'],'BA_output');

cell_size_mean = BA_output.whole_cell_size_mean;
total_vim_mean = BA_output.whole_cell_vim_totalamount_mean;
