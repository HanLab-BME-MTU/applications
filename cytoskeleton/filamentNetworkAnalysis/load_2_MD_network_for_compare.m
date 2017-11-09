function similarity_scoremap_cell = load_2_MD_network_for_compare(MD1_filename,MD2_filename,radius,show_save_everything_flag)
% function to compare two networks
% Input:   MD1_filename, MD2_filename:
%                       two MD file names,these two movie should be one channel with same number of frames.
%          radius:      the local neighborhood size for network comparison
%          show_save_everything_flag:
%                       whether to show and save all the figure during detailed process
%                       usually set to 0, to keep only the useful ones
% output:  similarity_scoremap_cell:
%          a cell structure with each frames similarity scoremap
%          every figure is saved to disc
%          at the end of the function the output dir is opened

% Liya Ding
% 2013

% input: show_save_everything_flag:
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

outdir = [MD_1.processes_{indexFilamentSegmentationProcess_1}.outFilePaths_{1},filesep,'similarity_results'];
if(~exist(outdir,'dir'))
    mkdir(outdir);
end


for iFrame = 1 : nFrame
    iFrame
    
    % load first movie
    MT_orientation = MD_1.processes_{indexFilamentSegmentationProcess_1}.loadChannelOutput(1,iFrame,'output','current_seg_orientation');
    MT_current_model = MD_1.processes_{indexFilamentSegmentationProcess_1}.loadChannelOutput(1,iFrame,'output','current_model');
    
    MT_current_seg = (isnan(MT_orientation)==0);
    
    MT_img =  MD_1.processes_{indexFlattenProcess_1}.loadChannelOutput(1,iFrame);
    %
    % load second movie
    VIF_orientation = MD_2.processes_{indexFilamentSegmentationProcess_2}.loadChannelOutput(1,iFrame+0,'output','current_seg_orientation');
    VIF_current_model = MD_2.processes_{indexFilamentSegmentationProcess_2}.loadChannelOutput(1,iFrame+0,'output','current_model');
    
    VIF_current_seg = (isnan(VIF_orientation)==0);
    
    VIF_img =  MD_2.processes_{indexFlattenProcess_2}.loadChannelOutput(1,iFrame);
    
    
    % % display the two channel together
    two_channel_img = zeros(size(VIF_img,1),size(VIF_img,2),3);
    two_channel_img(:,:,1) = VIF_img;
    two_channel_img(:,:,2) = MT_img;
    
    
    if(show_save_everything_flag==1)
        h1=figure(1);imagesc(two_channel_img/255);axis equal;axis off;
        saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.tif']);
        saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.fig']);
    end
    
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
    
    if(show_save_everything_flag==1)        
        % % plotting in white
        h2=figure(2);imagesc(two_channel_seg(:,1:end,:));axis equal;axis off;
        figure(2);
        set(gca, 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        
        saveas(h2,[outdir,filesep,'white_VIFMT_seg_frame_',num2str(iFrame),'.tif']);
        saveas(h2,[outdir,filesep,'white_VIFMT_seg_frame_',num2str(iFrame),'.fig']);
    end
     
    ch2_white_background_segment(find(VIF_current_seg==1))=0.85;
    two_channel_seg(:,:,2)= double(ch2_white_background_segment);
    
    if(show_save_everything_flag==1)        
        h2=figure(2);imagesc(two_channel_seg(:,1:end,:));axis equal;axis off;
        figure(2);
        saveas(h2,[outdir,filesep,'darkgreenVIF_white_VIFMT_seg_frame_',num2str(iFrame),'.tif']);
        saveas(h2,[outdir,filesep,'darkgreenVIF_white_VIFMT_seg_frame_',num2str(iFrame),'.fig']);
    end
     
    % getting the size of the images for comparison
    img_size = size(MT_img);
    
    % core function of the comparison
    [similarity_scoremap, difference_map]= network_similarity_scoremap(MT_current_model,VIF_current_model,img_size, radius);
    
    % save the output to disk
    
    save([outdir,filesep,'VIFMT_sm_maps_frame_',num2str(iFrame),'.mat'], ...
    'difference_map', 'similarity_scoremap');


    % plot the detailed results if requested.
    plot_differnce_map_wrapper(difference_map,outdir,iFrame,radius,show_save_everything_flag);
    
    % put results for all the frames together
    similarity_scoremap_cell{1, iFrame} = similarity_scoremap;
    
    %       close all;
end

winopen(outdir);
