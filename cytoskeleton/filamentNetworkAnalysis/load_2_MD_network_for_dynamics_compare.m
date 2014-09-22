function [similarity_scoremap_cell,similarity_scoremap_1to2_cell,similarity_scoremap_2to1_cell, distance_map_1_2_cell, distance_map_2_1_cell, angle_map_1_2_cell, angle_map_2_1_cell] ...
    = load_2_MD_network_for_dynamics_compare(MD1_filename,iChannel1, start_frame1, MD2_filename, iChannel2, start_frame2, radius,save_everything_flag)
% function to compare two networks
% Input:   MD1_filename,MD2_filename:
%                       two MD file names,these two movie should be one channel with same number of frames.
%          iChannel1, iChannel2:
%                       The channel in the movie to be compared.
%           start_frame1, start_frame2:
%                       the start frame of comparison; if 1 and 1, then
%                       same time; if 1 and 2 then the dynamics
%          radius:      the local neighborhood size for network comparison
%          save_everything_flag:
%                       whether to save all the figure during detailed process
%                       usually set to 0, to keep only the useful ones
% output:  similarity_scoremap_cell:
%          a cell structure with each frames similarity scoremap
%          every figure is saved to disc
%          at the end of the function the output dir is opened

% input: save_everything_flag:
MD1 = load(MD1_filename);
MD2 = load(MD2_filename);

indexMTChannel=1;
indexVIFChannel=1;

% check each of the MD files
MD = MD1.MD;
package_process_ind_script;
indexFilamentSegmentationProcess_1 = indexFilamentSegmentationProcess;
indexFlattenProcess_1 = indexFlattenProcess;
indexCellRefineProcess_1 = indexCellRefineProcess;
indexCellSegProcess_1 = indexCellSegProcess;

MD = MD2.MD;
package_process_ind_script;
indexFilamentSegmentationProcess_2 = indexFilamentSegmentationProcess;
indexFlattenProcess_2 = indexFlattenProcess;
indexCellRefineProcess_2 = indexCellRefineProcess;
indexCellSegProcess_2 = indexCellSegProcess;

MD_1 = MD1.MD;
MD_2 = MD2.MD;

nFrame = MD_1.nFrames_;

%initialize output
similarity_scoremap_cell = cell(1,nFrame);
similarity_scoremap_1to2_cell = cell(1,nFrame);
similarity_scoremap_2to1_cell = cell(1,nFrame);
      
distance_map_1_2_cell = cell(1,nFrame);
distance_map_2_1_cell = cell(1,nFrame);
angle_map_1_2_cell = cell(1,nFrame);
angle_map_2_1_cell = cell(1,nFrame);
         
outdir = [MD_1.processes_{indexFilamentSegmentationProcess_1}.outFilePaths_{iChannel1},filesep,'similarity_results'];
if(~exist(outdir,'dir'))
    mkdir(outdir);
end


for iFrame = start_frame1 : nFrame
    
    iFrame_2 = iFrame - start_frame1 + start_frame2;
    
    if iFrame_2<= nFrame
        iFrame
        
        % load first movie
        MT_orientation = MD_1.processes_{indexFilamentSegmentationProcess_1}.loadChannelOutput(iChannel1,iFrame,'output','current_seg_orientation');
        MT_current_model = MD_1.processes_{indexFilamentSegmentationProcess_1}.loadChannelOutput(iChannel1,iFrame,'output','current_model');
        
        MT_current_seg = (isnan(MT_orientation)==0);
        
        MT_img =  MD_1.processes_{indexFlattenProcess_1}.loadChannelOutput(iChannel1,iFrame);
        %
        % load second movie, at the defined frame number
        VIF_orientation = MD_2.processes_{indexFilamentSegmentationProcess_2}.loadChannelOutput(iChannel2,iFrame_2+0,'output','current_seg_orientation');
        VIF_current_model = MD_2.processes_{indexFilamentSegmentationProcess_2}.loadChannelOutput(iChannel2,iFrame_2+0,'output','current_model');
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        
        VIF_img =  MD_2.processes_{indexFlattenProcess_2}.loadChannelOutput(iChannel2,iFrame_2);
        
        % % display the two channel frame together
        two_channel_img = zeros(size(VIF_img,1),size(VIF_img,2),3);
        two_channel_img(:,:,1) = VIF_img;
        two_channel_img(:,:,2) = MT_img;
        
        h1=figure(1);imagesc(two_channel_img/255);axis equal;axis off;
        saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.tif']);
        saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.fig']);
        
        two_channel_seg= zeros(size(VIF_img,1),size(VIF_img,2),3);
        
        VIF_plus_MT_current_seg = VIF_current_seg + MT_current_seg >0;
        
        ch1_white_background_segment = double(MT_current_seg);
        ch2_white_background_segment = double(VIF_current_seg);
        ch3_white_background_segment = double(1-VIF_plus_MT_current_seg);
        
        ch1_white_background_segment(find(VIF_plus_MT_current_seg==0))=1;
        ch2_white_background_segment(find(VIF_plus_MT_current_seg==0))=1;
        
        two_channel_seg(:,:,1)= double(ch1_white_background_segment);
        two_channel_seg(:,:,2)= double(ch2_white_background_segment);
        two_channel_seg(:,:,3) = double(ch3_white_background_segment);
        
        % % plotting
        h2=figure(2);imagesc(two_channel_seg(:,1:end,:));axis equal;axis off;
        figure(2);
        set(gca, 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        
        %    saveas(h2,[outdir,filesep,'white_VIFMT_seg_frame_',num2str(iFrame),'.tif']);
        %    saveas(h2,[outdir,filesep,'white_VIFMT_seg_frame_',num2str(iFrame),'.fig']);
        
        ch2_white_background_segment(find(VIF_current_seg==1))=0.85;
        two_channel_seg(:,:,2)= double(ch2_white_background_segment);
        h2=figure(2);imagesc(two_channel_seg(:,1:end,:));axis equal;axis off;
        figure(2);
        saveas(h2,[outdir,filesep,'darkgreenVIF_white_VIFMT_seg_frame_',num2str(iFrame),'.tif']);
        saveas(h2,[outdir,filesep,'darkgreenVIF_white_VIFMT_seg_frame_',num2str(iFrame),'.fig']);
        
        img_size = size(MT_img);
        
        [similarity_scoremap, similarity_scoremap_1to2, similarity_scoremap_2to1,distance_map_1_2, distance_map_2_1, angle_map_1_2, angle_map_2_1] = network_similarity_scoremap(MT_current_model,VIF_current_model,img_size, radius,outdir,iFrame,save_everything_flag);
        
        [filament_similiarity_1,filament_similiarity_2] ...
            = filament_similarity(VIF_current_model,MT_current_model,img_size, ...
            radius,outdir,iFrame,save_everything_flag,...
            distance_map_1_2, distance_map_2_1, angle_map_1_2, angle_map_2_1)
        
        similarity_scoremap_cell{1, iFrame} = similarity_scoremap;
        similarity_scoremap_1to2_cell{1, iFrame} = similarity_scoremap_1to2;
        similarity_scoremap_2to1_cell{1, iFrame} = similarity_scoremap_2to1;
        distance_map_1_2_cell{1, iFrame} = distance_map_1_2;
        distance_map_2_1_cell{1, iFrame} = distance_map_2_1;
        angle_map_1_2_cell{1, iFrame} = angle_map_1_2;
        angle_map_2_1_cell{1, iFrame} = angle_map_2_1;
         
        %       close all;
    end
end


save([outdir,filesep,'VIFMT_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','similarity_scoremap_1to2_cell',...
    'similarity_scoremap_2to1_cell', ...
    'distance_map_1_2_cell','distance_map_2_1_cell', ...
    'angle_map_1_2_cell','angle_map_2_1_cell');


winopen(outdir);
