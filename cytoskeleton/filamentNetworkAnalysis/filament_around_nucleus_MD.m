function out_profiles = filament_around_nucleus_MD(MD, filament_channel, nucleus_channel, distance_pace,distance_max, display_flag)
% function to do filament density analysis wrt distance to nucleus with input MD, which has one
% channel as nucleus and one channel for filament

% input:    MD:                 the loaded movieData object for analysis
%           filament_channel:   the channel index containing filament
%           nucleus_channel:    the channel index containing nucleus
%           distance_pace:      the pace of calculating the layers as
%                               distance increases
%           distance_max:       the max distance to consider
%           display_flag:       flag, 1 if want to show; 0 if not to show

% output:   out_profiles is the output struct with these fields
%           int_mean_profile: the profiles telling when the increase of
%                   distance from the nucleus, what is the filament average intensity
%           st_mean_profile: the profiles telling when the increase of
%                   distance from the nucleus, what is the filament average
%                   steerable filter response
%           network_density_profile: the profiles telling when the increase of
%                   distance from the nucleus, what is the filament average
%                   density(as in segmented results)
%           int_median_profile: the  telling when the increase of
%                   distance from the nucleus, what is the filament
%                   intensity median
%           st_median_profile: the profile telling when the increase of
%                   distance from the nucleus, what is the filament
%                   steerable filter response median

%           Note all these above outputs, have the indices as
%           int_mean_profile(iFrame, iDistance)

%           the all frame outputs accordingly:

%           all_frame_int_mean_profile
%           all_frame_st_mean_profile
%           all_frame_network_mean_profile
%           all_frame_int_median_profile
%           all_frame_st_median_profile
%

% Liya Ding 2014

%%
% initialize the output varaiables, if no steerable filtering or filament
% segmentation, return empty for these output

out_profiles=[];

int_mean_profile=[];
st_mean_profile=[];
network_density_profile=[];
int_median_profile=[];
st_median_profile=[];

all_frame_int_mean_profile=[];
all_frame_st_mean_profile=[];
all_frame_network_density_profile=[];
all_frame_int_median_profile=[];
all_frame_st_median_profile=[];


%% Start the parameter/index gathering

movie_Dir = MD.outputDirectory_;
display_msg_flag=1;
package_process_ind_script;

% if no nucleus definition, return empty in everything
if(indexCellRefineProcess==0)
    display('need to do nucleus segmentation and mask refinement');
    return;
end

%% set up the output matrices

% the bins to examine
distance_bin = 0: distance_pace: distance_max;
num_distance_bin = numel(distance_bin)-1;

int_mean_profile = zeros(MD.nFrames_,num_distance_bin);
int_median_profile = zeros(MD.nFrames_,num_distance_bin);
all_frame_int_mean_profile = zeros(1,num_distance_bin);
all_frame_int_median_profile = zeros(1,num_distance_bin);
% pool through all frames
all_frame_int_pool = cell(1,num_distance_bin);

% if there is steerable filtering process
if(indexSteerabeleProcess>0)
    st_mean_profile = zeros(MD.nFrames_,num_distance_bin);
    st_median_profile = zeros(MD.nFrames_,num_distance_bin);
    all_frame_st_mean_profile = zeros(1,num_distance_bin);
    all_frame_st_median_profile = zeros(1,num_distance_bin);
    
    % pool through all frames
    all_frame_st_pool = cell(1,num_distance_bin);
    
end

% if there is filament segment process
if(indexFilamentSegmentationProcess>0)
    network_density_profile = zeros(MD.nFrames_,num_distance_bin);
    all_frame_network_density_profile = zeros(1,num_distance_bin);
    % pool through all frames
    all_frame_network_pool = cell(1,num_distance_bin);
end

for iD = 1 : num_distance_bin
    all_frame_int_pool{1,iD}=[];
    all_frame_st_pool{1,iD}=[];
    all_frame_network_pool{1,iD}=[];    
end

%% Start of the pool gathering the mean median calculation

for iFrame = 1:MD.nFrames_
    
    %% load image, cell segmentation, or steerable filtering and filament segmentation
    try
        filament_image = double(MD.channels_(filament_channel).loadImage(iFrame));
        nucleus_mask = MD.processes_{indexCellRefineProcess}.loadChannelOutput(nucleus_channel,iFrame);
    catch
        display('Cannot load the nucleus segmentation; or cannot find the filament image, pls check.');
        return;
    end
    
    % see if steerable filtering available
    if(indexSteerabeleProcess>0)        
        try
            st_image = MD.processes_{indexSteerabeleProcess}.loadChannelOutput(filament_channel,iFrame, 'output','MAX_st_res');
        catch
            st_image = [];
            display('Cannot load the steerable for the filament channel, pls check.');            
        end
    end
    
    
    % see if filament segmentation available
    if(indexFilamentSegmentationProcess>0)        
        try
            seg_image = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(filament_channel,iFrame, 'output','current_seg_orientation');
            seg_image = (isnan(current_seg_orientation)==0);
        catch
            seg_image=[];
            display('Cannot load the segmentation for the filament channel, pls check.');            
        end
    end
    
    
    dist_map = bwdist(nucleus_mask,'euclidean');
    
    for iD = 1 : num_distance_bin
        min_D = distance_bin(iD);
        max_D = distance_bin(iD+1);
        ind = find(dist_map>min_D & dist_map<=max_D);
        
        this_int_pool = filament_image(ind);
        int_mean_profile(iFrame, iD) = mean(this_int_pool);
        int_median_profile(iFrame, iD) = median(this_int_pool);
        
        all_frame_int_pool{1,iD} = [all_frame_int_pool{1,iD}; this_int_pool];
        
        % if not empty gather the st and seg
        if(~isempty(st_image))
            this_st_pool = st_image(ind);
            st_mean_profile(iFrame, iD) = mean(this_st_pool);
            st_median_profile(iFrame, iD) = median(this_st_pool);
            % all frame pool
            all_frame_st_pool{1,iD} = [all_frame_st_pool{1,iD}; this_st_pool];
        end
        
        if(~isempty(seg_image))
            this_seg_pool = seg_image(ind);
            network_density_profile(iFrame, iD) = mean(double(this_int_pool));
            % all frame pool
            all_frame_network_pool{1,iD} = [all_frame_network_pool{1,iD}; this_seg_pool];
        end
    end
    
    if(display_flag>0)
        figure(1);bar(int_mean_profile(iFrame,:));
        ylabel('Avg Int');
        xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
        title(['Frame ', num2str(iFrame),': Average intensity profiles as distance from nucleus increases']);
        axis auto;
        
        figure(2);bar(int_median_profile(iFrame,:));
        ylabel('Median Int');
        xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
        title(['Frame ', num2str(iFrame),': Median intensity profiles as distance from nucleus increases']);
        axis auto;
        
        
        if(~isempty(st_image))
            
            figure(3);bar(st_mean_profile(iFrame,:));
            ylabel('Avg St');
            xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
            title(['Frame ', num2str(iFrame),': Average steerable filtering result profiles as distance from nucleus increases']);
            axis auto;
            
            figure(4);bar(st_median_profile(iFrame,:));
            ylabel('Median st');
            xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
            title(['Frame ', num2str(iFrame),': Median steerable filtering result profiles as distance from nucleus increases']);
            axis auto;
        end
        
        if(~isempty(seg_image))
            figure(5);bar(network_density_profile(iFrame,:)*100);
            ylabel('Avg Filament Density, %');
            xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
            title(['Frame ', num2str(iFrame),': Average density profiles as distance from nucleus increases']);
            axis auto;
        end        
    end
end


for iD = 1 : num_distance_bin
    
    all_frame_int_mean_profile(1,iD) = mean(all_frame_int_pool{1,iD});
    all_frame_int_median_profile(1,iD) = median(all_frame_int_pool{1,iD});
    
    if(~isempty(all_frame_st_pool{1,iD}))        
        all_frame_st_mean_profile(1,iD) = mean(all_frame_st_pool{1,iD});
        all_frame_st_median_profile(1,iD) = median(all_frame_st_pool{1,iD});
    end
    
    if(~isempty(all_frame_network_pool{1,iD}))
        all_frame_network_density_profile(1,iD) = mean(double(all_frame_network_pool{1,iD}));
    end
end


 if(display_flag>0)
        figure(1);bar(all_frame_int_mean_profile);
        ylabel('Avg Int');
        xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
        title('All Frames: Average intensity profiles as distance from nucleus increases');
        axis auto;
        
        figure(2);bar(all_frame_int_median_profile);
        ylabel('Median Int');
        xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
        title('All Frames: Median intensity profiles as distance from nucleus increases');
        axis auto;
        
        
        if(~isempty(all_frame_st_mean_profile))
            
            figure(3);bar(all_frame_st_mean_profile);
            ylabel('Avg St');
            xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
            title('All Frames: Average steerable filtering result profiles as distance from nucleus increases');
            axis auto;
            
            figure(4);bar(all_frame_st_median_profile);
            ylabel('Median st');
            xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
            title('All Frames: Median steerable filtering result profiles as distance from nucleus increases');
            axis auto;
        end
        
        if(~isempty(all_frame_network_density_profile))
            figure(5);bar(all_frame_network_density_profile*100);
            ylabel('Avg Filament Density, %');
            xlabel(['Distance from nuclues , in each ', num2str(distance_pace),' pixel bins']);
            title('All Frames: Average density profiles as distance from nucleus increases');
            axis auto;
        end
        
    end

%% gather outputs
out_profiles=[];

out_profiles.int_mean_profile = int_mean_profile;
out_profiles.st_mean_profile =  st_mean_profile;
out_profiles.network_density_profile =  network_density_profile;
out_profiles.int_median_profile = int_median_profile;
out_profiles.st_median_profile = st_median_profile;

out_profiles.all_frame_int_mean_profile = all_frame_int_mean_profile;
out_profiles.all_frame_st_mean_profile =  all_frame_st_mean_profile;
out_profiles.all_frame_network_density_profile =all_frame_network_density_profile;
out_profiles.all_frame_int_median_profile = all_frame_int_median_profile;
out_profiles.all_frame_st_median_profile = all_frame_st_median_profile;


