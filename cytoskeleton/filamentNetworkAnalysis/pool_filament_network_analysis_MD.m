function network_feature_pool_cell = pool_filament_network_analysis_MD(MD)
% this function is for pooling the network analysis results for all frame
% in a movie
% Input:  MD: loaded MD with network analysis done already
% Output:  network_feature_pool_cell:  a cell structure in the length of the
%                        channel numbers, each have a structure that are
%                        pooled results from all frames in this channel

%          the fields are:
%                    straightness_per_filament_pool_all 
%                    orientation_pixel_pool_display_all 
%                    length_per_filament_pool_all 
%                    pixel_number_per_filament_pool_all 
%                    density_filament_all 
%                    scrabled_density_filament_all 
%                    intensity_per_filament_pool_all 
%                    mean_intensity_per_filament_pool_all 
%                    intensity_per_fat_filament_pool_all 
%                    mean_intensity_per_fat_filament_pool_all 
%                    scale_per_filament_pool_all 
%     
% Created 07 2014 by Liya Ding, Matlab R2012b

movie_Dir = MD.outputDirectory_;

% find all the index of different processes
package_process_ind_script;

% initialize the cell structures
network_feature = cell(length(MD.channels_),nFrame);
network_feature_pool_cell = cell(length(MD.channels_),1);

%% loading
% for each channel do analysis
for iChannel = 1 :  length(MD.channels_)
    
    % define where the results are stored
    outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
    
    % for each frame
    for iFrame = 1 : nFrame
        
        % load from mat file
        load([outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'network_analysis.mat'], ...
            'output_feature');
        
        % put them into the cell for future pooling
        network_feature{iChannel,iFrame} = output_feature;
    end
end

%% Pooling
% for each channel do analysis
for iChannel = 1 :  length(MD.channels_)
    
    network_feature_pool_cell{iChannel,1} = [];
    network_feature_pool_cell{iChannel,1}.straightness_per_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.orientation_pixel_pool_display_all = [];
    network_feature_pool_cell{iChannel,1}.length_per_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.pixel_number_per_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.density_filament_all = [];
    network_feature_pool_cell{iChannel,1}.scrabled_density_filament_all = [];
    network_feature_pool_cell{iChannel,1}.intensity_per_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.mean_intensity_per_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.intensity_per_fat_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.mean_intensity_per_fat_filament_pool_all = [];
    network_feature_pool_cell{iChannel,1}.scale_per_filament_pool_all = [];
    
    % for each frame
    for iFrame = 1 : nFrame
        
        output_feature =  network_feature{iChannel,iFrame};
        
        network_feature_pool_cell{iChannel,1}.straightness_per_filament_pool_all = [ network_feature_pool_cell{iChannel,1}.straightness_per_filament_pool_all; ...
            output_feature.straightness_per_filament_pool(:)];
        network_feature_pool_cell{iChannel,1}.orientation_pixel_pool_display_all = [network_feature_pool_cell{iChannel,1}.orientation_pixel_pool_display_all;...
            output_feature.orientation_pixel_pool_display(:)];
        
        network_feature_pool_cell{iChannel,1}.length_per_filament_pool_all = [network_feature_pool_cell{iChannel,1}.length_per_filament_pool_all;...
            output_feature.length_per_filament_pool(:) ];
        network_feature_pool_cell{iChannel,1}.pixel_number_per_filament_pool_all = [network_feature_pool_cell{iChannel,1}.pixel_number_per_filament_pool_all ;...
            output_feature.pixel_number_per_filament_pool(:) ];
        network_feature_pool_cell{iChannel,1}.density_filament_all = [network_feature_pool_cell{iChannel,1}.density_filament_all;...
            output_feature.density_filament ];
        network_feature_pool_cell{iChannel,1}.scrabled_density_filament_all = [network_feature_pool_cell{iChannel,1}.scrabled_density_filament_all;...
            output_feature.scrabled_density_filament ];
        network_feature_pool_cell{iChannel,1}.intensity_per_filament_pool_all = [network_feature_pool_cell{iChannel,1}.intensity_per_filament_pool_all;...
            output_feature.intensity_per_filament_pool(:) ];
        network_feature_pool_cell{iChannel,1}.mean_intensity_per_filament_pool_all = [network_feature_pool_cell{iChannel,1}.mean_intensity_per_filament_pool_all;...
            output_feature.mean_intensity_per_filament_pool(:) ];
        network_feature_pool_cell{iChannel,1}.intensity_per_fat_filament_pool_all = [network_feature_pool_cell{iChannel,1}.intensity_per_fat_filament_pool_all;...
            output_feature.intensity_per_fat_filament_pool(:) ];
        network_feature_pool_cell{iChannel,1}.mean_intensity_per_fat_filament_pool_all = [network_feature_pool_cell{iChannel,1}.mean_intensity_per_fat_filament_pool_all;...
            output_feature.mean_intensity_per_fat_filament_pool(:) ];
        network_feature_pool_cell{iChannel,1}.scale_per_filament_pool_all = [network_feature_pool_cell{iChannel,1}.scale_per_filament_pool_all;...
            output_feature.scale_per_filament_pool(:) ];
        
    end
    
    %     % make all of them col vectors(except the density maps)
    %     network_feature_pool_cell{iChannel,1}.straightness_per_filament_pool_all = network_feature_pool_cell{iChannel,1}.straightness_per_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.orientation_pixel_pool_display_all = network_feature_pool_cell{iChannel,1}.orientation_pixel_pool_display_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.length_per_filament_pool_all = network_feature_pool_cell{iChannel,1}.length_per_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.pixel_number_per_filament_pool_all = network_feature_pool_cell{iChannel,1}.pixel_number_per_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.intensity_per_filament_pool_all = network_feature_pool_cell{iChannel,1}.intensity_per_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.mean_intensity_per_filament_pool_all = network_feature_pool_cell{iChannel,1}.mean_intensity_per_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.intensity_per_fat_filament_pool_all = network_feature_pool_cell{iChannel,1}.intensity_per_fat_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.mean_intensity_per_fat_filament_pool_all = network_feature_pool_cell{iChannel,1}.mean_intensity_per_fat_filament_pool_all(:);
    %
    %     network_feature_pool_cell{iChannel,1}.scale_per_filament_pool_all = network_feature_pool_cell{iChannel,1}.scale_per_filament_pool_all(:);
    %
end
