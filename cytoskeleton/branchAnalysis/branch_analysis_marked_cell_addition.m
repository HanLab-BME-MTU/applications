function [BA_output, f_display] = branch_analysis_marked_cell_addition(MD, iChannel, iCell, half_size,min_branch_size_Threshold,filament_stat_flag,figure_flag)
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

if(sum(isCompleted)<10)    
    return;
end

if(~exist(outputPath,'dir'))
    mkdir(outputPath);
end

CompletedFrame = find(isCompleted>0);
CompletedFrame = CompletedFrame(:);
CompletedFrame = CompletedFrame';

% in this addition version, load existing result

% keep the frame numbers
BA_output.CompletedFrame = CompletedFrame;

% these index is sure to be continuous
FirstFrame = min(CompletedFrame);
nCompleteFrame = length(CompletedFrame);

raw_mask_cell = cell(1,nCompleteFrame);
smoothed_mask_cell = cell(1,nCompleteFrame);
current_seg_cell= cell(1,nCompleteFrame);
orienation_map_filtered_cell= cell(1,nCompleteFrame);
region_orientation_cell = cell(1,nCompleteFrame);

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
    % default as empty in filament segmentation
    current_seg=[];
    orienation_map_filtered=[];

    % if filament stat is available and requested 
    if(exist(GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'steerable_vote_',filename_short_strs{iFrame},'.mat']),'file') && filament_stat_flag>0)
        load( GetFullPath([ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(VIF_channel),filesep,'DataOutput',filesep,'steerable_vote_',...
        filename_short_strs{iFrame},'.mat']),'current_seg','orienation_map_filtered');
    end
     
     current_seg_cell{1,iCompleteFrame} = current_seg;
     orienation_map_filtered_cell{1,iCompleteFrame} = orienation_map_filtered;

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
    
    current_mask = raw_mask_cell{1,iCompleteFrame};
    
    current_seg = current_seg_cell{1,iCompleteFrame};
    
    smoothed_current_mask = smoothed_mask_cell{1,iCompleteFrame};
    
    try
        C_xy = regionprops(smoothed_current_mask,'Centroid');
        center_x(iCompleteFrame) = C_xy.Centroid(1);
        center_y(iCompleteFrame) = C_xy.Centroid(2);
    catch
        disp('error in cell mask centroid');
        break;
    end
    
    RG_framem1 = zeros(size(smoothed_current_mask,1),size(smoothed_current_mask,2),3)>0;
    RG_framem1(:,:,1) = (smoothed_current_mask);
    if iCompleteFrame>1
        RG_framem1(:,:,2) = (smoothed_mask_cell{1,iCompleteFrame-1});
    end
    
    boundary_map = bwmorph(smoothed_current_mask,'remove');

   thin_current_mask = bwmorph(smoothed_current_mask,'thin',Inf);
    
    % Find the end points and branching points
    end_poins_map = bwmorph(thin_current_mask,'endpoint');
    
    branching_points_map = bwmorph(thin_current_mask,'branchpoints');
    end_points_map = bwmorph(thin_current_mask,'endpoints');
    branching_points_map = imdilate(branching_points_map,ones(3,3));
    % Delete these branching points for now
    
    skel_no_branching = (thin_current_mask - branching_points_map)>0;
     
     branching_points_map = bwmorph(skel_no_branching,'branchpoints');
     branching_points_map = imdilate(branching_points_map,ones(3,3));
    
    skel_no_branching = (skel_no_branching - branching_points_map)>0;
    
    % not % break the cener as well
    
    % Label all isolated lines(curves)
    [labelMask, nL] = bwlabel(skel_no_branching);
    
    % Get properties for each of curve
    ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid','Orientation');
    
    
    % Redefine variable for easy of notation
    obAreas = [ob_prop.Area];
    
    
     obOrientation = ([ob_prop.Orientation]).*pi/180;
    
    region_orientation = nan(size(labelMask));
    
    nLine = length(obAreas);
    
    % Compute dist map
    [D,IDX] = bwdist(skel_no_branching);
    
    region_branch_label = zeros(size(skel_no_branching));
    branch_leaf_flag=[];
    region_branch_wo_inner_label = region_branch_label;
    
    labelMask_without_inner = labelMask;
    for iL = 1 : nL
        this_branch_IDX = find(D==0 & labelMask==iL);
        % if there is no ending point, then it is a fake branch
        
        Lia = ismember(IDX(:),this_branch_IDX);
       
         if(sum(end_points_map(find(labelMask==iL)))==0 && ...
                sum(boundary_map(find(Lia>0)))==0 && ...
                length(this_branch_IDX)<min_branch_size_Threshold)
            branch_leaf_flag(iL)=0;
            labelMask_without_inner(find(labelMask==iL))=0;
            region_branch_wo_inner_label(Lia>0)=0;
        else
            branch_leaf_flag(iL)=1;   
            region_branch_label(Lia>0)=iL;   
            region_branch_wo_inner_label(Lia>0)=iL;              
         end
        
         
         region_orientation(Lia>0)=obOrientation(iL);
        region_branch_label_R(Lia>0)=(0.3+color_array(iL,1))/(1.3);
        region_branch_label_G(Lia>0)=(0.3+color_array(iL,2))/(1.3);
        region_branch_label_B(Lia>0)=(0.3+color_array(iL,3))/(1.3);
    end
    
    region_branch_label(smoothed_current_mask==0)=0;
%     region_branch_wo_inner_label(smoothed_current_mask==0)=0;
   
    region_branch_label_cell{iCompleteFrame} = uint8(region_branch_label);
%     region_branch_wo_inner_label_cell{iCompleteFrame} = region_branch_wo_inner_label;
    
    label_skel_cell{iCompleteFrame} = uint8(labelMask);
%     label_skel_wo_inner_cell{iCompleteFrame} = labelMask_without_inner;
    
    branch_leaf_flag_cell{iCompleteFrame} = branch_leaf_flag;
    branch_number_per_frame(iCompleteFrame)=sum(branch_leaf_flag);
    
    region_orientation_cell{iCompleteFrame} =region_orientation;
end


movieInfo =[];

BA_output.cell_travel_length = ...
    sum(sqrt((center_x(2:end)-center_x(1:end-1)).^2 + (center_y(2:end)-center_y(1:end-1)).^2 ));

BA_output.cell_travel_distance = ...
    (sqrt((center_x(end)-center_x(1)).^2 + (center_y(end)-center_y(1)).^2 ));

BA_output.cell_travel_speed = ...
    BA_output.cell_travel_length./nCompleteFrame;

BA_output.cell_marked_frame_number = ...
    nCompleteFrame;

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

h11=figure(11);hold off;
plot(center_x,center_y,'r');
hold on;
plot(smooth_center_x,smooth_center_y);
% take the atan2, negative y due to the y axis is top small bottom big in
% image matrix
trajectory_angle = atan2(-(smooth_center_y(2:end)-smooth_center_y(1:end-1)),smooth_center_x(2:end)-smooth_center_x(1:end-1));

saveas(h11,[outputPath,filesep,'path_smoothed.tif']);


%%
for iCompleteFrame = 1 : nCompleteFrame
    
    labelMask = label_skel_cell{iCompleteFrame};
    
    branch_leaf_flag = branch_leaf_flag_cell{iCompleteFrame};

    % get rid of inner fake branches
    for iB = 1 : length(branch_leaf_flag)
        if branch_leaf_flag(iB)==0;
            labelMask(labelMask==iB)=0;
        end        
    end
    
    % Get properties for each of curve
    ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid','Orientation');
    
    % Redefine variable for easy of notation
    obAreas = [ob_prop.Area];
    obAreas = obAreas';
    
    obCentroid = zeros(2, length(obAreas));
    obCentroid(:) = [ob_prop.Centroid];
    obCentroid = obCentroid';
    
    obOrientation = ([ob_prop.Orientation]).*pi/180;
    
%     region_orientation = nan(size(labelMask));
%     
%     for iR = 1 : length(obAreas)
%     
%         region_orientation(labelMask==iR) = obOrientation(iR);
%     
%     end
    
    movieInfo(iCompleteFrame).xCoord = [obCentroid(:,1) obCentroid(:,1)*0];
    
    movieInfo(iCompleteFrame).yCoord = [obCentroid(:,2) obCentroid(:,2)*0];
    
    movieInfo(iCompleteFrame).zCoord = [obAreas*0 obAreas*0];
   
%     movieInfo(iCompleteFrame).Orientation = obOrientation;
    
    movieInfo(iCompleteFrame).amp = [obAreas obAreas*0];
       
%     region_orientation_cell{1, iCompleteFrame} = region_orientation;
    
end

gapCloseParam.mergeSplit=0;
gapCloseParam.minTrackLen=1;
gapCloseParam.timeWindow = 0;

costMatrices(1).funcName = 'costMatLinearMotionLink';

parameters.linearMotion = 0; %no linear motion
parameters.minSearchRadius = 20.5;%1.5; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = 100.0;%1.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = 1; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

costMatrices(1).parameters = parameters;


costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';

%parameters

%needed all the time
parameters.linearMotion = 0;%1; %use linear motion Kalman filter.

parameters.minSearchRadius = 1;%1.5; %minimum allowed search radius.
parameters.maxSearchRadius = 2.0;%1.5; %maximum allowed search radius.
parameters.brownStdMult = 3*1; %multiplication factor to calculate Brownian search radius from standard deviation.
parameters.timeReachConfB = 2; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).

parameters.ampRatioLimit = [0.5 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = 1; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 3*ones(1,1); %multiplication factor to calculate linear search radius from standard deviation.
parameters.timeReachConfL = 1; %same as timeReachConfB, but for the linear part of the motion.
parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

costMatrices(2).parameters = parameters;

kalmanFunctions.reserveMem = 'kalmanResMemLM';
kalmanFunctions.initialize = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

try
[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(...
    movieInfo,costMatrices,gapCloseParam,kalmanFunctions,[],...
    [],[]);
catch    
    disp('Error in tracking');
    return;
end

new_label_skel_cell =  label_skel_cell;
new_region_branch_label_cell = region_branch_label_cell;
for iCompleteFrame = 1 : nCompleteFrame
    new_label_skel_cell{iCompleteFrame} = new_label_skel_cell{iCompleteFrame}*0;
    new_region_branch_label_cell{iCompleteFrame} = new_region_branch_label_cell{iCompleteFrame}*0;    
end

% final output of the number of branches tracked
BA_output.branch_number_tracked = length(tracksFinal);

% the maximum number branches among all frames
BA_output.branch_number_max = max(branch_number_per_frame);

% the mean number branches among all frames
BA_output.branch_number_mean = mean(branch_number_per_frame);


for iT = 1 : length(tracksFinal)
    track_ind = tracksFinal(iT).tracksFeatIndxCG;
    track_start = tracksFinal(iT).seqOfEvents(1,1);
    track_end = tracksFinal(iT).seqOfEvents(2,1);
    
    for iF = track_start :track_end
        this_skel = label_skel_cell{iF};
        new_skel = new_label_skel_cell{iF};
        new_skel(this_skel==(track_ind(iF-track_start+1)))=iT;
        new_label_skel_cell{iF} = new_skel;
        
        this_region = region_branch_label_cell{iF};
        new_region =  new_region_branch_label_cell{iF};
        new_region(this_region==(track_ind(iF-track_start+1)))=iT;
        new_region_branch_label_cell{iF} = new_region;
   end
    BA_output.branch_duration_array(iT) = track_end - track_start+1;
end
BA_output.branch_duration_mean = mean(BA_output.branch_duration_array);


color_array = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 1 0; 0 1 1; rand(max(494,max(length(tracksFinal))-6),3)];
      

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more analysis and plotting
branch_analysis_after_tracking;

%%
% display and save results

f_display = figure();

set(f_display,'menu','none','toolbar','none');
fid = fopen([outputPath,filesep,'branch_analysis_report.txt'],'wt');
display_message_window = uipanel(f_display,'Units','normalized', ...
    'position',[0.1 0.1 0.8 0.8], ...
    'title',['Branch Analysis for Channel ', num2str(iChannel), ', Cell ', num2str(iCell),':'],...
    'FontSize',12);

str_line = ['In the ', num2str(BA_output.cell_marked_frame_number),' marked frames, '];
fprintf(fid, [str_line,'\n  \r\n']);
strings{1}=str_line;

str_line = ['this cell has its center travelled ', num2str(BA_output.cell_travel_length),' pixels.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{2}=str_line;

str_line = ['The cell center moved at speed ', num2str(BA_output.cell_travel_speed),' pixels/frame.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{3}=str_line;

str_line = ['The cell center moved ', num2str(BA_output.cell_travel_distance),' pixels.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{4}=str_line;

str_line = ['The branches tracked: ', num2str(BA_output.branch_number_tracked),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{5}=str_line;

str_line = ['The maximum number branches among all frames: ', ...
    num2str(BA_output.branch_number_max),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{6}=str_line;


str_line = ['The average number branches: ', ...
    num2str(BA_output.branch_number_mean),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{7}=str_line;

str_line = ['The average duration of branches: ', ...
    num2str(BA_output.branch_duration_mean),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{8}=str_line;


str_line = ['The average VIM intensity in protrusion region: ', ...
    num2str(BA_output.protrusion_vif_mean_intensity),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{9}=str_line;

str_line = ['The average VIM intensity in retraction region: ', ...
    num2str(BA_output.retraction_vif_mean_intensity),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{10}=str_line;

str_line = ['The average VIM intensity in the whole cell: ', ...
    num2str(BA_output.whole_cell_vif_mean_intensity),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{11}=str_line;

str_line = ['The average VIM Total Amount in the whole cell: ', ...
    num2str(BA_output.whole_cell_vim_totalamount_mean),'.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{12}=str_line;

str_line = ['The average cell size: ', ...
    num2str(BA_output.whole_cell_size_mean),' pixels.'];
fprintf(fid, [str_line,'\n  \r\n']);
strings{13}=str_line;

fclose(fid);
 
lbh = uicontrol(display_message_window,'style','text','Units','normalized','position',...
    [0 0 1 1],'FontSize',12);
set(lbh,'HorizontalAlignment','left');

set(lbh,'string',strings);

save([outputPath,filesep,'branch_analysis_results.mat'],'BA_output');

