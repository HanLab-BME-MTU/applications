function network_feature = load_MD_network_for_analysis(MD,ROI,radius,figure_flag, save_everything_flag,feature_flag)
% function to do network analysis with input MD

% input:    MD:    the loaded movieData object.
%           ROI:   the user input ROI, if not defined [], will use the whole
%                  image area
%           radius: the definition of neighborhood
%           figure_flag: 1 to plot histgrams, 0 not to
%           save_everything_flag: 1 to save the plot, 0 not to
%           feature_flag: the flag for controlling which feature to
%                         calculate, a vector of 16 bits, 1 for need to
%                         calculate, 0 for skip it to save time.

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

% output:   network_feature, a cell structure for each channel, each frame.
%           Each struct with field of the 16 features as above
            
% Liya Ding 2013

% if no input, ROI is full image
if(nargin<2)
    ROI = [];
end

% if no input for radius, set it as default 20
if(nargin<3)
    radius = 20;
end

% if no input to plot or not, plot
if(nargin<4)
    figure_flag = 1;
end

% if no input to save plot or not, save
if(nargin<5)
    save_everything_flag = 1;
end

% if no input as which feature to calculate, do for all features
if(nargin<6)
    feature_flag = ones(1,16);
end
%% Get movie data ready

movie_Dir = MD.outputDirectory_;

% find all the index of different processes
display_msg_flag = 0; % display warning or not
package_process_ind_script;
network_feature=cell(length(MD.channels_),nFrame);

% for each channel do analysis
for iChannel = 1 :  length(MD.channels_)
    
    outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
     
    SteerableChannelOutputDir = MD.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    
    % make out put directory if not existing
    if(~exist(outdir,'dir'))
        mkdir(outdir);
    end
    Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % for each frame
    for iFrame = 1 : nFrame
        display(['iChannel: ', num2str(iChannel),', iFrame:', num2str(iFrame)]);
        %% % Load the data
       
        % load scale map
         % this line in commandation for shortest version of filename
        filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
            
        try
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_short_strs{iFrame},'.mat']);            
        catch
            % in the case of only having the short-old version
            if nFrame>1
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_shortshort_strs{iFrame},'.mat']);            
            else
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_shortshort_strs,'.mat']);            
            end
        end
        
        % load the filament segmentation results       
        VIF_tip_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','tip_orientation');
        VIF_RGB_seg_orient_heat_map = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','RGB_seg_orient_heat_map');
        
        VIF_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
        VIF_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
        
%         %% % this part is for tip analysis, not ready, just ignore this part
%         % find the tip location
%         [vim_tip_y,vim_tip_x] = find(~isnan(VIF_tip_orientation));
%         
%         % load the cell segmentation results from segmentation and mask
%         % refine processes
%         CellSegmentation = MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame);
%         
%         % define the cell boundary
%         CellBoundary = bwboundaries(CellSegmentation);
%         CellBoundary = CellBoundary{1}; %% x in second col
%         
%         % display the image with the tips
%         h1 = figure(1);
%         imagesc(VIF_RGB_seg_orient_heat_map);axis off; axis image;
%         hold on;
%         plot(CellBoundary(:,2),CellBoundary(:,1),'r.');
%         plot(vim_tip_x,vim_tip_y,'*');
        
        %% % transfer into digital representation and define ROI
        
        [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
            = filament_model_to_digital_with_orientation(VIF_current_model);
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        
        % if the input ROI is [], then use the whole area
        if(isempty(ROI))
            ROI = ones(size(VIF_current_seg));
        end
 
        min_length = MD.processes_{indexFilamentSegmentationProcess}.funParams_.LengthThreshold;

        if(numel(min_length)>1)
        min_length = min_length(iChannel);
        end
        
        %% % do analysis
        im_name = MD.channels_(iChannel).getImageFileNames(iFrame);
        current_img =  MD.channels_(iChannel).loadImage(iFrame);
                
        % get network feature that is only related to the network
        [output_network_features, VIF_ROI_model, VIF_ROI_orientation_model] ...
            = network_analysis(VIF_current_model,...
                VIF_current_seg, ROI, radius,feature_flag);
        
        % get network feature that is related to intensity or st image
         output_image_features = perfilament_analysis(VIF_ROI_model, VIF_ROI_orientation_model,...
                 VIF_current_seg, scaleMap,current_img, MAX_st_res, feature_flag);
        
        % putting the network features together
        output_feature =  cat_struct(output_network_features,output_image_features);
        
        % plot the network features in hists
        network_features_plotting(output_feature, figure_flag, save_everything_flag, feature_flag,...
                im_name, outdir,iChannel,iFrame)
%         close all;

        % add one last component, the cell_mask
        if(indexCellRefineProcess>0)
            Cell_Mask = ROI.*((MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame))>0);
        else
            Cell_Mask = ROI;
        end
        output_feature.Cell_Mask = Cell_Mask;
        
        % save output feature for single image(single channel, single frame)
        save([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.mat'],...
            'output_feature');
     
        % put output feature to cell for all channels, all frames
        network_feature{iChannel,iFrame} = output_feature;
        
    end
    
    % save output feature for all channels, all frames)       
    save([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_allframe.mat'],...
            'network_feature');
end

 
        