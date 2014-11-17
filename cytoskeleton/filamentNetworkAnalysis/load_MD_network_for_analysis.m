function network_feature = load_MD_network_for_analysis(MD,ROI,radius)
% function to do network analysis with input MD

% input:    MD:    the loaded movieData object.
%           ROI:   the user input ROI, if not defined [], will use the whole
%                  image area
% output:   network_feature, a cell structure for each channel, each frame.
%           Each struct with field of:
%           'straightness_per_filament_pool','orientation_pixel_pool_display',...
%           'length_per_filament_pool'. 
            
% Liya Ding 2013

movie_Dir = MD.outputDirectory_;

% find all the index of different processes
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
                
        output_feature = network_analysis(VIF_current_model,VIF_orientation, ...
            VIF_current_seg,outdir,iChannel,iFrame,ROI,min_length,...
            im_name,radius);
        
        output_feature_Filaments = perfilament_analysis(VIF_current_model,VIF_orientation, ...
            VIF_current_seg,outdir,iChannel,iFrame,ROI,min_length,...
            im_name,radius,scaleMap,current_img,MAX_st_res);
        
        output_feature=  cat_struct(output_feature,output_feature_Filaments);
                    
%         close all;

        Cell_Mask = ROI.*((MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame))>0);
          
        output_feature.Cell_Mask = Cell_Mask;
            save([outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.mat'],...
            'output_feature');
     
        network_feature{iChannel,iFrame} =output_feature;
        
    end
    save([outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_allframe.mat'],...
            'network_feature');
end

 
        