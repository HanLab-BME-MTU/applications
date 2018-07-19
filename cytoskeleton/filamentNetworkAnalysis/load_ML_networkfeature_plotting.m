function load_ML_networkfeature_plotting(ML,mode_index, figure_flag, save_everything_flag, feature_flag)
% function to plot network analysis results with input MD

% input:    ML:    the loaded movieList object.
%           mode_index:   the mode for analysis:
%                         1 : to tell what is the pooled distribution for
%                         the whole movie list
%                         2 : to tell which one within this list is
%                         different from the others
%
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

% the number of movies
nMovie =  length(ML.movieDataFile_);
networkanalysis_output_ML_cell= cell(1,1);

MD = MovieData.load(ML.movieDataFile_{1});

% all the movie should have the same number of channels, for same content
% each channel
nChannel = length(MD.channels_);

% in order to reduce memory need,  work channel by channel

for iChannel = 1 :  nChannel
    
     clearvars -except 'nMovie' 'iChannel' 'nChannel' 'networkanalysis_output_ML_cell' 'ML' 'figure_flag' 'save_everything_flag' 'feature_flag'
        
     close all;
    
     %%  initialize all possible pools
        
    if mode_index == 1
        straightness_per_filament_pool_all = [];
        length_per_filament_pool_all = [];
        pixel_number_per_filament_pool_all = [];
        density_filament_all = [];
        scrabled_density_filament_all = [];
        orientation_pixel_pool_display_all = [];
        orientation_pixel_pool_display_center_all = [];
        intensity_per_filament_pool_all = [];
        mean_intensity_per_filament_pool_all = [];
        intensity_per_fat_filament_pool_all = [];
        mean_intensity_per_fat_filament_pool_all = [];
        scale_per_filament_pool_all = [];
        st_per_filament_pool_all = [];
        mean_st_per_filament_pool_all =[];
        st_per_fat_filament_pool_all =[];
        mean_st_per_fat_filament_pooll_all =[];
    else
        % for mode 2, only the mean of the each distribtion is taken
        quarters_straightness_per_filament_pool_all = [];
        quarters_length_per_filament_pool_all = [];
        quarters_pixel_number_per_filament_pool_all = [];
        quarters_density_filament_all = [];
        quarters_scrabled_density_filament_all = [];
        quarters_orientation_pixel_pool_display_all = [];
        quarters_orientation_pixel_pool_display_center_all = [];
        quarters_intensity_per_filament_pool_all = [];
        quarters_mean_intensity_per_filament_pool_all = [];
        quarters_intensity_per_fat_filament_pool_all = [];
        quarters_mean_intensity_per_fat_filament_pool_all = [];
        quarters_scale_per_filament_pool_all = [];
        quarters_st_per_filament_pool_all = [];
        quarters_mean_st_per_filament_pool_all =[];
        quarters_st_per_fat_filament_pool_all =[];
        quarters_mean_st_per_fat_filament_pooll_all =[];        
    end
    
    %%  run through all movies for current channel
    
    for iM  = 1 : nMovie
              
        % load this movie
        MD = MovieData.load(ML.movieDataFile_{iM});
        % the number of channels
        display('======================================================================');
        display(['iM:', num2str(iM)]);
        
        %% Get movie data ready
        
        movie_Dir = MD.outputDirectory_;
        
        % find all the index of different processes
        display_msg_flag = 0; % display warning or not
        package_process_ind_script;
                
        outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
        
        % for each frame
        for iFrame = 1 : nFrame
            display(['iChannel: ', num2str(iChannel),', iFrame:', num2str(iFrame)]);
            %% % Load the network analysis results
            load([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.mat'],...
                'output_feature');
            
            im_name = MD.channels_(iChannel).getImageFileNames(iFrame);
            
            if mode_index == 1
                %% in mode 1, telling whole distribution
                %  so pool every requested feature
                
                % in order to prevent matrix to oversized, creat a
                % subsample rate, rate 1 for feature, rate 2 for pixel
                % based distributions
                
                % assume length is always calculated....
                subsample_index_1 = randsample(length(output_feature.length_per_filament_pool), ...
                    length(output_feature.length_per_filament_pool)/(max(1, nMovie/10)));
                
                %  subsample_index_2 only if pixel based features are
                %  calculated
                %                 subsample_index_2 = randsample(length(output_feature.density_filament), ...
                %                     length(output_feature.density_filament)/(max(1, nMovie/5)));
                
                % output according to user request
                
                if(feature_flag(1)>0)
                    straightness_per_filament_pool_all = [straightness_per_filament_pool_all;
                        (output_feature.straightness_per_filament_pool(subsample_index_1))'];
                end
                
                if(feature_flag(2)>0)
                    length_per_filament_pool_all = [length_per_filament_pool_all;
                        (output_feature.length_per_filament_pool(subsample_index_1))'];
                end
                
                if(feature_flag(3)>0)
                    pixel_number_per_filament_pool_all = [ pixel_number_per_filament_pool_all;
                        (output_feature.pixel_number_per_filament_pool(subsample_index_1))'];
                end
                
                
                if(feature_flag(4)>0 ||  feature_flag(5)>0 )
                    % get masked density information
                    density_masked = output_feature.density_filament(output_feature.Cell_Mask>0);
                    
                    subsample_index_2 = randsample(numel(density_masked), ...
                        numel(density_masked)/(max(1, nMovie/5)));
                    
                    density_filament_all = [ density_filament_all;
                        density_masked(subsample_index_2)];
                    
                    try
                        scrabled_density_masked = output_feature.scrabled_density_filament(output_feature.Cell_Mask>0);
                        
                        scrabled_density_filament_all = [ scrabled_density_filament_all;
                            scrabled_density_masked(subsample_index_2)];
                    end
                end
                
                if(feature_flag(6)>0 || feature_flag(7)>0)
                    
                    subsample_index_3 = randsample(length(output_feature.orientation_pixel_pool_display), ...
                        length(output_feature.orientation_pixel_pool_display)/(max(1, nMovie/8)));
                    
                    orientation_pixel_pool_display_all = [ orientation_pixel_pool_display_all;
                        output_feature.orientation_pixel_pool_display(subsample_index_3)];
                    
                    try
                        orientation_pixel_pool_display_center_all = [ orientation_pixel_pool_display_center_all;
                            output_feature.orientation_pixel_pool_display_center(subsample_index_2)];
                    end
                end
                
                if(feature_flag(8)>0)
                    intensity_per_filament_pool_all = [intensity_per_filament_pool_all;
                        output_feature.intensity_per_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(9)>0)
                    mean_intensity_per_filament_pool_all = [mean_intensity_per_filament_pool_all;
                        output_feature.mean_intensity_per_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(10)>0)
                    intensity_per_fat_filament_pool_all = [intensity_per_fat_filament_pool_all;
                        output_feature.intensity_per_fat_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(11)>0)
                    mean_intensity_per_fat_filament_pool_all = [mean_intensity_per_fat_filament_pool_all;
                        output_feature.mean_intensity_per_fat_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(12)>0)
                    scale_per_filament_pool_all = [scale_per_filament_pool_all;
                        output_feature.scale_per_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(13)>0)
                    st_per_filament_pool_all = [st_per_filament_pool_all;
                        output_feature.st_per_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(14)>0)
                    mean_st_per_filament_pool_all = [mean_st_per_filament_pool_all;
                        output_feature.mean_st_per_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(15)>0)
                    st_per_fat_filament_pool_all = [st_per_fat_filament_pool_all;
                        output_feature.st_per_fat_filament_pool(subsample_index_1)];
                end
                
                if(feature_flag(16)>0)
                    mean_st_per_fat_filament_pool_all = [mean_st_per_fat_filament_pool_all;
                        output_feature.mean_st_per_fat_filament_pool(subsample_index_1)];
                end
                
            else
                %% in the difference telling mode
                if mode_index == 2
                    % with this script, get the 25%,50%,and 75% for each
                    % distribution
                    quarter_pooling_network_features;
                end
            end
        end
        
    end
end


