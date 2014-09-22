function twochannel_network_comparison(MD)
% this function is to compare the filaments in two channels, in distance
% and orientation
% Input: the the movieData object

movie_Dir = MD.getPath;

% find the information of this MD and the index of processed ran
package_process_ind_script;

% assume comparing VIF(channel1) and MT(channel2)
indexVIFChannel = 1;
indexMTChannel = 2;  

outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{1},filesep,'similarity_results'];

if(~exist(outdir,'dir'))
    mkdir(outdir);
end

dist_pool_for_crossing=[];
ang_pool_for_crossing=[];
    
for iFrame = 13 : 29
    iFrame
    
    % VIF_current_model
    % load from the segmentation
    VIF_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(2,iFrame+0,'output','current_model');

    VIF_img =  MD.processes_{indexFlattenProcess}.loadChannelOutput(2,iFrame);
   
    % MT models
    MT_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(1,iFrame,'output','current_model');
    % load the image
    MT_img =  MD.processes_{indexFlattenProcess}.loadChannelOutput(1,iFrame);
    

    % image size
    img_size = size(VIF_img);
   
    % get the bw segmentation
    VIF_current_seg = filament_model_to_seg_bwim(VIF_current_model,img_size,[]);
    MT_current_seg = filament_model_to_seg_bwim(MT_current_model,img_size,[]);

   
    % two channel images
    two_channel_img = zeros(size(VIF_img,1),size(VIF_img,2),3);
    two_channel_img(:,:,1)=VIF_img;
    two_channel_img(:,:,2)=MT_img;
    h1=figure(1);imagesc(two_channel_img/255);axis equal;axis off;
    saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.tif']);
    
    % two channel segmentation
    two_channel_seg= zeros(size(VIF_img,1),size(VIF_img,2),3);
    two_channel_seg(:,:,1)=VIF_current_seg;
    two_channel_seg(:,:,2)=MT_current_seg;
    
    h2=figure(2);imagesc(two_channel_seg);axis equal;axis off;
    saveas(h2,[outdir,filesep,'VIFMT_seg_frame_',num2str(iFrame),'.tif']);
    
    radius=60;
    
    similarity_scoremap = network_similarity_scoremap(VIF_current_model,MT_current_model,img_size, radius,outdir,iFrame);
    
end
