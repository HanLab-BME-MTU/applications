movie_Dir = '/home/ld94/files/LCCB/interfil/Liya/fromJess/2013/130522_VimMT/130522_VimMT_WTandNoc_Analyze_Drive/130522_VimMT_RPE_60x_002/crop_1';

load([movie_Dir, filesep, 'movieData.mat']);

nFrame = MD.nFrames_;

VIF_model = cell(1,nFrame);
MT_model = cell(1,nFrame);
flatten_dir{1} = MD.processes_{3}.outFilePaths_{1};
flatten_dir{2} = MD.processes_{3}.outFilePaths_{2};
outdir = [MD.processes_{5}.outFilePaths_{1},filesep,'similarity_results'];
    mkdir(outdir);
    
    dist_pool_for_crossing=[];
    ang_pool_for_crossing=[];
    
for iFrame = 1 : 10
    iFrame
    VIF_orientation = MD.processes_{5}.loadChannelOutput(1,iFrame+0,'output','current_seg_orientation');
    VIF_current_model = MD.processes_{5}.loadChannelOutput(1,iFrame+0,'output','current_model');
    
%     VIF_current_model
    [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
    = filament_model_to_digital_with_orientation(VIF_current_model);

%     VIF_orientation = VIF_orientation(140:200,150:250);
   
    VIF_current_seg = (isnan(VIF_orientation)==0);
    VIF_img =  imread([flatten_dir{1},filesep,'flatten_',num2str(iFrame+0,'%03d'),'.tif']);
    
    MT_orientation = MD.processes_{5}.loadChannelOutput(2,iFrame,'output','current_seg_orientation');
%     MT_orientation = MT_orientation(140:200,150:250);
    MT_current_model = MD.processes_{5}.loadChannelOutput(2,iFrame,'output','current_model');
    
    [MT_digital_model,MT_orientation_model,MT_XX,MT_YY,MT_OO] ...
        = filament_model_to_digital_with_orientation(MT_current_model);

    
    MT_current_seg = (isnan(MT_orientation)==0);
    
    MT_img =  imread([flatten_dir{2},filesep,'flatten_',num2str(iFrame,'%03d'),'.tif']);
%     MT_img = MT_img(140:200,150:250);
%    VIF_img = VIF_img(140:200,150:250);
   
    
    
    two_channel_img = zeros(size(VIF_img,1),size(VIF_img,2),3);
    two_channel_img(:,:,1)=VIF_img;
    two_channel_img(:,:,2)=MT_img;
    h1=figure(1);imagesc(two_channel_img/255);axis equal;axis off;
     saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.tif']);
  saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.fig']);
   
    
    
    two_channel_seg= zeros(size(VIF_img,1),size(VIF_img,2),3);
    two_channel_seg(:,:,1)=VIF_current_seg;
     two_channel_seg(:,:,2)=MT_current_seg;
    
    h2=figure(2);imagesc(two_channel_seg);axis equal;axis off;
   saveas(h2,[outdir,filesep,'VIFMT_seg_frame_',num2str(iFrame),'.tif']);
   saveas(h2,[outdir,filesep,'VIFMT_seg_frame_',num2str(iFrame),'.fig']);
   
    network_similarity_scoremap;
    
      close all;
end
