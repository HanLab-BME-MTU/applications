function [network_feature,network_feature_MD_allCh_wholepool_cell] = load_MD_network_for_analysis(MD,CellROI,radius,figure_flag, save_everything_flag,feature_flag,vimscreen_flag)
% function to do network analysis with input MD

% input:    MD:    the loaded movieData object.
%           CellROI:   the user input CellROI, if not defined [], will use the whole
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

%           17:   output_feature.filament_mean_curvature
%           18:   output_feature.curvature_per_pixel_pool 

%           19:   output_feature.profileCell
%           20:   output_feature.profileAllCell

%           21:   output_feature.Centripetal_fila 
%           22:   output_feature.Centripetal_pixel 


% 
% textcell = {'Straightness of filament (for each filament)',... % 1
%     'Length (for each filament)',...                           % 2
%     'Pixel number of segmentation (for each filament)',...     % 3
%     'Filament density (for each valid pixel)',...                   % 4
%     'Scrabled filament density (for each valid pixel)',...          % 5
%     'Filament orientation (for each pixel)',...                % 6
%     'Filament orientation (for each pixel,centered)',...       % 7
%     'Filamet intensity (integrated for each filament)',...           % 8
%     'Filamet intensity (average for each filament)',...              % 9
%     'Filamet intensity (integrated for each dilated filament)',...   % 10
%     'Filamet intensity (average for each dilated filament)',...      % 11
%     'Scale of detected filaments (for each pixel)',...                   % 12
%     'ST Responce (integrated for each filament)',...            % 13
%     'ST Responce (average for each filament)',...               % 14
%     'ST Responce (integrated for each dilated filament)',...    % 15
%     'ST Responce (average for each dilated filament)',...       % 16
%     'Filament curvature (average for each filament)',...     % 17
%     'Filament curvature (for each pixel on filaments)',...   % 18
%     'Filament Profiles for each cell (for vim screen only)',... %19   
%     'Filament Profiles for all cells (for vim screen only)',... %20
%     'Filament Centripetal angle for each filament (for vim screen only)',... %21 
%     'Filament Centripetal angle for each pixel on filament (for vim screen only)',... %22
%     'Cell Mask',... % 23
%     'Number of Nucleus'... %24
%     }


% output:   network_feature, a cell structure for each channel, each frame.
%           Each struct with field of the 16 features as above
            
% Liya Ding 2013

% if no input, CellROI is full image
if(nargin<2)
    CellROI = [];
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
    feature_flag = ones(1,18);
end

if(numel(feature_flag)<18)
    feature_flag = [feature_flag(:); zeros(18-numel(feature_flag),1)];
end

% if no input as if it is vim screen, set to no
if(nargin<7)
    vimscreen_flag = 0;
end


%% Get movie data ready

movie_Dir = MD.outputDirectory_;
wholemovie_output_dir = [MD.outputDirectory_,filesep,'whole_movie_network_analysis'];
 
if(~exist(wholemovie_output_dir,'dir'))
    mkdir(wholemovie_output_dir);
end

% find all the index of different processes
display_msg_flag = 0; % display warning or not
package_process_ind_script;
network_feature=cell(length(MD.channels_),nFrame);

if(vimscreen_flag>0 && (indexFilamentSegmentationProcess==0 ||indexFlattenProcess==0 ...
        ||indexSteerabeleProcess==0 || indexFilamentPackage==0) )
    MD = vimscreen_forceMDAddProcessSave(MD);
    package_process_ind_script;
end


% for each channel that did filament segmentation, do analysis
try
    validChannels = MD.processes_{indexFilamentSegmentationProcess}.funParams_.ChannelIndex;
catch
    validChannels = 1 : numel(MD.channels_);
end

% into row vector
validChannels = validChannels(:);
validChannels = validChannels'; % into row vector

for iChannel = validChannels
    
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
        frame_tic = tic;
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
        
        %% % transfer into digital representation and define CellROI
        
        [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
            = filament_model_to_digital_with_orientation(VIF_current_model);
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        
        % if the input CellROI is [], then use the whole area
        if(isempty(CellROI))
            CellROI = ones(size(VIF_current_seg));
        end
                
        % see if there is cell segmentation to rule out noisy background
        Cell_Mask = CellROI;
        try
            Cell_Mask = CellROI.*((MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame))>0);                       
        end
        
 
        min_length = MD.processes_{indexFilamentSegmentationProcess}.funParams_.LengthThreshold;

        if(numel(min_length)>1)
        min_length = min_length(iChannel);
        end
        
        %% % do analysis
        im_name = MD.channels_(iChannel).getImageFileNames(iFrame);
        current_img =  MD.channels_(iChannel).loadImage(iFrame);
                
        % get network feature that is only related to the network
        display(' --- network features');
        tic
        [output_network_features, VIF_ROI_model, VIF_ROI_orientation_model] ...
            = network_analysis(VIF_current_model,...
            VIF_current_seg, Cell_Mask, radius,feature_flag);
        toc
        
        display(' --- intensity, scale, steerable-response features');   
        tic
        % get network feature that is related to intensity or st image
        output_image_features = perfilament_analysis(VIF_ROI_model, VIF_ROI_orientation_model,...
            VIF_current_seg, scaleMap,current_img, MAX_st_res, feature_flag);
        toc
        
        % putting the network features together
        output_feature =  cat_struct(output_network_features,output_image_features);
        
        
        %% for vim screen the third channel is nucleus, count the number of cells
        if(vimscreen_flag>0 && indexCellRefineProcess >0)
           display(' --- Vim Screen - Structure Features');   
        
            % load the nucleus channel segmentation
            nucleus_mask = (MD.processes_{indexCellRefineProcess}.loadChannelOutput(3,iFrame))>0;
            
            % label the nucleus segmentation image and get the number of
            % cells(nucleus)
           [labelMaskNucleus, numNucleus] = bwlabel(nucleus_mask);
           output_feature.number_of_nucleus = numNucleus;
           
           % With the nucleus locations, get the vim properties in each
           % devided region of each cell
           tic
           vim_output_feature = vim_screen_network_features(labelMaskNucleus,...
               VIF_current_seg,VIF_ROI_model,current_img, nms,feature_flag,40);           
           toc
           % put these into the final struct
           output_feature =  cat_struct(output_feature,vim_output_feature);        
        end           
        
        % add one last component, the cell_mask
        Cell_Mask = CellROI;
        if(indexCellRefineProcess>0 && vimscreen_flag == 0)
            try
                Cell_Mask = CellROI.*((MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame))>0);
              Cell_Mask(Cell_Mask==0)=nan;
      
            end            
        end
        
        
        output_feature.Cell_Mask = Cell_Mask;
        
        if(~isempty(Cell_Mask))
            
            output_feature.density_filament =  (output_feature.density_filament).*(output_feature.Cell_Mask);
            output_feature.scrabled_density_filament =  (output_feature.scrabled_density_filament).*(output_feature.Cell_Mask);
        end
        
        output_feature.filament_density_mean = nanmean(output_feature.density_filament(:));
        output_feature.scrabled_density_filament_mean = nanmean(output_feature.scrabled_density_filament(:));
        
         
%         % save output feature for single image(single channel, single frame)
%         save([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.mat'],...
%             'output_feature');
%         
        % new name save output feature for single image(single channel, single frame)
        save([outdir,filesep,'network_feature_ch_',num2str(iChannel),'_',filename_short_strs{iFrame},'.mat'],...
            'output_feature');
        
        
        tic
        % plot the network features in hists
        network_features_plotting(output_feature, figure_flag, save_everything_flag, feature_flag,vimscreen_flag,...
                im_name, outdir,iChannel,iFrame)
%         close all;
        toc
        
        % put output feature to cell for all channels, all frames
        network_feature{iChannel,iFrame} = output_feature;
        
       display(' this frame totaling: ');        
       toc(frame_tic)
    end
    
    network_feature_MD_thisCh_wholepool = network_feature_pool_MD_gather(network_feature, iChannel, feature_flag);
    network_feature_MD_allCh_wholepool_cell{iChannel} = network_feature_MD_thisCh_wholepool;
    
    network_features_plotting(network_feature_MD_allCh_wholepool_cell, figure_flag, save_everything_flag, feature_flag,vimscreen_flag,...
        im_name, wholemovie_output_dir, iChannel);
        
    % save output feature for all channels(till this), all frames)       
    save([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_allframe.mat'],...
            'network_feature','network_feature_MD_thisCh_wholepool');
end
   
% save output feature for all channels, all frames)
save([wholemovie_output_dir,filesep,'network_analysis_feature_allch_allframe.mat'],...
    'network_feature','network_feature_MD_allCh_wholepool_cell');
