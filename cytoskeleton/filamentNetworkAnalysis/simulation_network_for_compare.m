ROOT_DIR = 'C:\Ding\MT_data\Break_data';
groundtruth_mat_dir = [ROOT_DIR, filesep, 'MAT'];

movie_Dir = [ROOT_DIR, filesep, 'B', filesep, '0.5'];

load([movie_Dir, filesep, 'movieData.mat']);

indexVIFChannel = 1;
indexMTChannel = 1;

package_process_ind_script;

SM_model = cell(1,nFrame);
MT_model = cell(1,nFrame);

outdir = [MD.outputDirectory_,filesep,'evaluation'];
mkdir(outdir);
    
    dist_pool_for_crossing=[];
    ang_pool_for_crossing=[];
    
for iFrame = 1 : 1
    iFrame
    
    im_filename = MD.channels_(1).getImageFileNames(iFrame);
    im_filename = im_filename{1};
    
    % get the ID of this simulated image
    ID_filename = im_filename(1:end-4);
    
    % load corresponding mat file for the ground truth and break lists
    load([groundtruth_mat_dir, filesep, ID_filename,'.mat']);
    
    
%     afterrelease_MT_body_x_cell
%     afterrelease_MT_body_y_cell
%     break_list
%     
    SM_current_model = cell(length(afterrelease_MT_body_x_cell),1);
    
    for iMT = 1 : length(afterrelease_MT_body_x_cell)
        SM_current_model{iMT,1}(:,1) = (afterrelease_MT_body_x_cell{iMT})';
        SM_current_model{iMT,1}(:,2) = (afterrelease_MT_body_y_cell{iMT})';
    end
    MT_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(1,iFrame,'output','current_model');

    MT_img = MD.channels_(1).loadImage(iFrame);
    
    img_size = size(MT_img);
    
    % get the bw segmentation
    SM_current_seg = filament_model_to_seg_bwim(SM_current_model,img_size,[]);
    MT_current_seg = filament_model_to_seg_bwim(MT_current_model,img_size,[]);

    two_channel_seg= zeros(size(MT_img,1),size(MT_img,2),3);
    two_channel_seg(:,:,1)=SM_current_seg;
    two_channel_seg(:,:,2)=MT_current_seg;
    
    h2=figure(2);imagesc(two_channel_seg);axis equal;axis off;
    saveas(h2,[outdir,filesep,'SM_MT_seg_frame_',num2str(iFrame),'.tif']);
   
    network_similarity_with_breaks;
    
    close all;
end
