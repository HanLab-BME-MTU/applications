function load_MD_networkfeature_plotting(MD,figure_flag, save_everything_flag,feature_flag)
% function to plot network analysis results with input MD

% input:    MD:    the loaded movieData object.
%           figure_flag: 1 to plot histgrams, 0 not to
%           save_everything_flag: 1 to save the plot, 0 not to
%           feature_flag: the flag for controlling which feature to
%                         plot, a vector of 16 bits, 1 for need to
%                         plot, 0 for skip it to save time.

%           %% here is a list of feature_flag corrspondence:
%
%            1:   output_feature.straightness_per_filament_pool
%            2:   output_feature.length_per_filament_pool
%            3:   output_feature.pixel_number_per_filament_pool
%            4:   output_feature.density_filament
%            5:   output_feature.scrabled_density_filament
%
%            6:   output_feature.orientation_pixel_pool_display
%            7:   output_feature.orientation_pixel_pool_display_center

%            8:   output_feature.intensity_per_filament_pool
%            9:   output_feature.mean_intensity_per_filament_pool
%           10:   output_feature.intensity_per_fat_filament_pool
%           11:   output_feature.mean_intensity_per_fat_filament_pool
%           12:   output_feature.scale_per_filament_pool
%
%           13:   output_feature.st_per_filament_pool
%           14:   output_feature.mean_st_per_filament_pool
%           15:   output_feature.st_per_fat_filament_pool
%           16:   output_feature.mean_st_per_fat_filament_pool

% output:   plots of network_feature

% Liya Ding 2014

% if no input to plot or not, plot
if(nargin<2)
    figure_flag = 1;
end

% if no input to save plot or not, save
if(nargin<3)
    save_everything_flag = 1;
end

% if no input as which feature to calculate, do for all features
if(nargin<4)
    feature_flag = ones(1,16);
end
%% Get movie data ready

movie_Dir = MD.outputDirectory_;

% find all the index of different processes
display_msg_flag = 0; % display warning or not
package_process_ind_script;

% for each channel do analysis
for iChannel = 1 :  length(MD.channels_)
    
    outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
    
    % for each frame
    for iFrame = 1 : nFrame
        display(['iChannel: ', num2str(iChannel),', iFrame:', num2str(iFrame)]);
        %% % Load the network analysis results
        load([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.mat'],...
            'output_feature');
        
        im_name = MD.channels_(iChannel).getImageFileNames(iFrame);
        
        % plot the network features in hists
        network_features_plotting(output_feature, figure_flag, save_everything_flag, feature_flag,...
            im_name, outdir,iChannel,iFrame)
        %         close all;
    end
end

