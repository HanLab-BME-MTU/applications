function network_feature =...
    load_MD_network_for_VIF_anisotropy(MD,CellROI,radius,...
    figure_flag, save_everything_flag,feature_flag,set_visible)
% function to do network analysis with input MD

% input:    MD:    the loaded movieData object.
%           CellROI:   the user input CellROI, if not defined [], will use the whole
%                  image area
%           radius: the definition of neighborhood
%           figure_flag: 1 to plot histgrams, 0 not to
%           save_everything_flag: 1 to save the plot, 0 not to
%           feature_flag: the flag for controlling which feature to
%                         calculate, a vector of 2, 1 for need to
%                         calculate, 0 for skip it to save time.

%           %% here is a list of feature_flag corrspondence:
%
%            1:   output_feature.filament_orientation_STD
%            2:   output_feature.filament_orientation_STD_median
%
% textcell = {'Pool of Filament orientation STD',... % 1
%     'Median of Pool of Filament orientation STD ',...     % 2
%     }

% output:   network_feature, a cell structure for each channel, each frame.
%           Each struct with field of the 2 features as above

% Liya Ding 2015

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

% if no input as which feature to calculate, do for all first 16 features
if(nargin<6)
    feature_flag = ones(1,2);
end

% if no input as if display figure or not, display them
if(nargin<7)
    set_visible = 1;
end

if(~exist('save_tif_flag','var'))
    save_tif_flag = 0;
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
network_feature_MD_allCh_wholepool_cell=cell(1,length(MD.channels_));

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
    
%     try
%         outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
%     catch
%         continue;
%     end
%     
%     % make out put directory if not existing
%     if(~exist(outdir,'dir'))
%         mkdir(outdir);
%     end
    Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % for each frame
    for iFrame = 1 : nFrame
        display(['iChannel: ', num2str(iChannel),', iFrame:', num2str(iFrame)]);
        frame_tic = tic;
        %% % Load the data        
        try
            % load the filament segmentation results
            VIF_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
            VIF_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
        catch
            % if there is no segmentation
            continue;
        end
        
        % load scale map
        % this line in commandation for shortest version of filename
        filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        VIF_current_seg_dilated = imfilter(VIF_current_seg, ones(radius*2+1,radius*2+1));
        
        % if the input CellROI is [], then use the whole area
        if(isempty(CellROI))
            CellROI = ones(size(VIF_current_seg));
        end
        
        % see if there is cell segmentation to rule out noisy background
        Cell_Mask = CellROI;
        try
            Cell_Mask = CellROI.*((MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame))>0);
        end
        
        Cell_Mask = Cell_Mask.*VIF_current_seg_dilated;
        
        Cell_Mask(1:radius,:)=0;
        Cell_Mask(:,1:radius)=0;
        Cell_Mask(end-radius:end,:)=0;
        Cell_Mask(:,end-radius:end)=0;        
        
         %% % transfer into digital representation and define CellROI
        if(exist('VIF_current_model','var'))
          [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
                = filament_model_to_digital_with_orientation(VIF_current_model);
            VIF_orientation_modified = nan(size(VIF_orientation));
            VIF_orientation_modified(sub2ind(size(VIF_orientation),VIF_YY,VIF_XX)) = VIF_OO;
            VIF_orientation = VIF_orientation_modified;
        end
        
         %% % do analysis
        % get network feature that is only related to the network
        display(' --- network features');
        
        filament_orientation_STD_Map = nan(size(VIF_orientation));
        
        [indy,indx]=find(Cell_Mask>0);
        
        for indi = 1 : length(indx)
            y=indy(indi);
            x=indx(indi);
            patch_O = VIF_orientation(y-radius:y+radius,x-radius:x+radius);
            patch_O(patch_O>pi/2) = patch_O(patch_O>pi/2) - pi;
            patch_O(patch_O>pi/2) = patch_O(patch_O>pi/2) - pi;
            patch_O(patch_O<=-pi/2) = patch_O(patch_O<=-pi/2) + pi;
            patch_O(patch_O<=-pi/2) = patch_O(patch_O<=-pi/2) + pi;
            
            if(numel(patch_O(isnan(patch_O)==0))>3)  
                if(x==205&&y==607)
                    stophere=1;
                end
                filament_orientation_STD_Map(y,x) = abs(circ_std(patch_O(isnan(patch_O)==0)*2)/2);            
            end
        end
        
        output_feature = [];
        output_feature.filament_orientation_STD = filament_orientation_STD_Map(isnan(filament_orientation_STD_Map)==0);
        output_feature.filament_orientation_STD_median = median(output_feature.filament_orientation_STD);
        
        % new name save output feature for single image(single channel, single frame)
        save([wholemovie_output_dir,filesep,'network_feature_vifstd_ch_',num2str(iChannel),'_',filename_short_strs{iFrame},'.mat'],...
            'output_feature');
                
        % put output feature to cell for all channels, all frames
        network_feature{iChannel,iFrame} = output_feature;
        
        display(' this frame totaling: ');
        toc(frame_tic)
    end 
end

try
    % save output feature for all channels, all frames)
    save([wholemovie_output_dir,filesep,'network_analysis_feature_vifstd_allch_allframe.mat'],...
        'network_feature');
end
