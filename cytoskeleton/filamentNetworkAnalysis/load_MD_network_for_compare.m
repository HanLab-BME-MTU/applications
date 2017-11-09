function [similarity_scoremap_cell,difference_map_cell]=...
    load_MD_network_for_compare(MD,radius,show_save_everything_flag,...
    longest_radius,sigma_gaussian, sigma_d, sigma_theta)
% function to compare two networks
% Input:   MD: movieData object loaded
%                       this movie should be two channel movie
%          radius:      the local neighborhood size for network comparison
%          save_everything_flag:
%                       whether to save all the figure during detailed process
%                       usually set to 0, to keep only the useful ones
%          longest_radius: the longest distance to be considered, here put a big number so that later process could be done
%          sigma_d, sigma_theta: if not the default set here
% output:  similarity_scoremap_cell:
%          a cell structure with each frames similarity scoremap
%          every figure is saved to disc
%          at the end of the function the output dir is opened


% if no input for longest radius, set it as default 100
if(nargin<4)
    longest_radius = 100;
end

flag_default = 0;

if(nargin<5)
    sigma_gaussian = 3*radius/8;
    flag_default = 1;
end
if(nargin<6)
    sigma_d = sqrt(3)*radius/4;
end

if(nargin<7)
    sigma_theta = pi/(2*sqrt(3));
end

if(isempty(sigma_gaussian))
     sigma_gaussian = 3*radius/8;
end

if(isempty(sigma_d))
     sigma_d = sqrt(3)*radius/4;
end

if(isempty(sigma_theta))
     sigma_theta = pi/(2*sqrt(3));
end

sigma_gaussian_ratio =     sigma_gaussian/((3*radius/8));
sigma_d_ratio =     sigma_d/((sqrt(3)*radius/4)) ;
sigma_theta_ratio =     sigma_theta/((pi/(2*sqrt(3))));
 

package_process_ind_script;

movie_Dir = MD.outputDirectory_;

nFrame = MD.nFrames_;

VIF_model = cell(1,nFrame);
MT_model = cell(1,nFrame);
% flatten_dir{1} = MD.processes_{3}.outFilePaths_{1};
% flatten_dir{2} = MD.processes_{3}.outFilePaths_{2};
outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),'_similarity_results'];

if(~exist(outdir,'dir'))
    mkdir(outdir);
end

if(flag_default==1)
    customized_outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),...
        '_similarity_r',num2str(radius),'_s_default'];
else
      customized_outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),...
        '_similarity_r',num2str(radius),'_sg',num2str(round(sigma_gaussian_ratio*100),'%d'),...
        '_sd',num2str(round(sigma_d_ratio*100),'%d'),...
        '_sa',num2str(round(sigma_theta_ratio*100),'%d')];
%     customized_outdir(customized_outdir=='.')='p';
end


if(~exist(outdir,'dir'))
    mkdir(outdir);
end
if(~exist(customized_outdir,'dir'))
    mkdir(customized_outdir);
end

dist_pool_for_crossing=[];
ang_pool_for_crossing=[];


%initialize output
similarity_scoremap_cell = cell(1,nFrame);
difference_map_cell = cell(1,nFrame);

intensity_pool_VIF = [];
intensity_pool_MT = [];

for iFrame = 1 : nFrame
    iFrame
    
    MT_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(1,iFrame,'output','current_seg_orientation');
    MT_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(1,iFrame,'output','current_model');
    
    [MT_digital_model,MT_orientation_model,MT_XX,MT_YY,MT_OO] ...
        = filament_model_to_digital_with_orientation(MT_current_model);
    
    MT_current_seg = (isnan(MT_orientation)==0);
    
    %     MT_img =  MD.processes_{indexFlattenProcess}.loadChannelOutput(1,iFrame);
    MT_img =  MD.channels_(1).loadImage(iFrame);
    
    VIF_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(2,iFrame+0,'output','current_seg_orientation');
    VIF_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(2,iFrame+0,'output','current_model');
    
    %     VIF_current_model
    [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
        = filament_model_to_digital_with_orientation(VIF_current_model);
    
    VIF_current_seg = (isnan(VIF_orientation)==0);
    
    %   VIF_img =  MD.processes_{indexFlattenProcess}.loadChannelOutput(2,iFrame);
    VIF_img =  MD.channels_(2).loadImage(iFrame);
    
    
    % show the two image segmentation together
    two_channel_img = zeros(size(VIF_img,1),size(VIF_img,2),3);
    two_channel_img(:,:,1) = MT_img;
    two_channel_img(:,:,2) = VIF_img;
    
    if(show_save_everything_flag==1)
        h1=figure(1);imagesc(double(two_channel_img)./double(max(max(max(MT_img)), max(max(VIF_img)))));
        axis equal;axis off;
        saveas(h1,[outdir,filesep,'MTVIF_img_frame_',num2str(iFrame),'.tif']);
        saveas(h1,[outdir,filesep,'MTVIF_img_frame_',num2str(iFrame),'.fig']);
    end
    %
    two_channel_seg= zeros(size(VIF_img,1),size(VIF_img,2),3);
    
    
    VIF_plus_MT_current_seg = VIF_current_seg + MT_current_seg >0;
    
    intensity_pool_VIF = [intensity_pool_VIF VIF_img(VIF_plus_MT_current_seg>0)];
    intensity_pool_MT = [intensity_pool_MT MT_img(VIF_plus_MT_current_seg>0)];
    
    ch1_white_background_segment = double(MT_current_seg);
    ch2_white_background_segment = double(VIF_current_seg);
    ch3_white_background_segment = double(1-VIF_plus_MT_current_seg);
    
    ch1_white_background_segment(find(VIF_plus_MT_current_seg==0))=1;
    ch2_white_background_segment(find(VIF_plus_MT_current_seg==0))=1;
    
    two_channel_seg(:,:,1)= double(ch1_white_background_segment);
    two_channel_seg(:,:,2)= double(ch2_white_background_segment);
    two_channel_seg(:,:,3) = double(ch3_white_background_segment);
    
    if(show_save_everything_flag==1)
        
        h2=figure(2);imagesc(two_channel_seg(:,1:end,:));axis equal;axis off;
        figure(2);
        %     set(gca, 'Position', get(gca, 'OuterPosition') - ...
        %        get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        
        saveas(h2,[outdir,filesep,'white_MTVIF_seg_frame_',num2str(iFrame),'.tif']);
        saveas(h2,[outdir,filesep,'white_MTVIF_seg_frame_',num2str(iFrame),'.fig']);
    end
    
    ch2_white_background_segment(find(VIF_current_seg==1))=0.85;
    
    two_channel_seg(:,:,2)= double(ch2_white_background_segment);
    if(show_save_everything_flag==1)
        
        h2=figure(2);imagesc(two_channel_seg(:,1:end,:));axis equal;axis off;
        figure(2);
        saveas(h2,[outdir,filesep,'darkgreenVIF_white_MTVIF_seg_frame_',num2str(iFrame),'.tif']);
        saveas(h2,[outdir,filesep,'darkgreenVIF_white_MTVIF_seg_frame_',num2str(iFrame),'.fig']);
    end
    
    img_size = size(MT_img);
    
    % core function of the comparison
    [similarity_scoremap, difference_map]= network_similarity_scoremap(MT_current_model,VIF_current_model,img_size, radius,...
        longest_radius,sigma_d, sigma_theta,sigma_gaussian);
    
    % save the output to disk
    
    save([outdir,filesep,'VIFMT_sm_maps_frame_',num2str(iFrame),'.mat'], ...
        'difference_map', 'similarity_scoremap');
   save([customized_outdir,filesep,'VIFMT_sm_maps_frame_',num2str(iFrame),'.mat'], ...
        'difference_map', 'similarity_scoremap');
    
    % plot the detailed results if requested.
    plot_differnce_map_wrapper(difference_map,outdir,customized_outdir,iFrame,radius,show_save_everything_flag);
    
    % put results for all the frames together
    similarity_scoremap_cell{1, iFrame} = similarity_scoremap;
    difference_map_cell{1, iFrame} = difference_map;
    
    
    save([outdir,filesep,'Similarity_maps_frame',num2str(iFrame),'.mat'], 'similarity_scoremap', 'difference_map');
    save([customized_outdir,filesep,'Similarity_maps_frame',num2str(iFrame),'.mat'], 'similarity_scoremap', 'difference_map');
    
    
    %       close all;
end

        
save([outdir,filesep,'VIFMT_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','difference_map_cell');
        
save([customized_outdir,filesep,'VIFMT_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','difference_map_cell',...
    
% 
% winopen(outdir);
% 
