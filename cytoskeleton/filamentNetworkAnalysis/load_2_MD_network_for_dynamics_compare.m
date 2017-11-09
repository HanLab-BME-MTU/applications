function [similarity_scoremap_cell,similarity_scoremap_1to2_cell,similarity_scoremap_2to1_cell, distance_map_1_2_cell, distance_map_2_1_cell, angle_map_1_2_cell, angle_map_2_1_cell] ...
    = load_2_MD_network_for_dynamics_compare(MD1_filename,iChannel1, start_frame1,...
    MD2_filename, iChannel2, start_frame2, ...
    radius,show_save_everything_flag,...
    longest_radius,sigma_gaussian, sigma_d, sigma_theta,save_tif_flag)
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

if(nargin<9)
    longest_radius = radius;
end

flag_default = 0;

if(nargin<10)
    sigma_gaussian = 3*radius/8;
    flag_default = 1;
end
if(nargin<11)
    sigma_d = sqrt(3)*radius/4/2;
end

if(nargin<12)
    sigma_theta = pi/(2*sqrt(3))/2;
end

if(isempty(sigma_gaussian))
     sigma_gaussian = 3*radius/8;
end

if(isempty(sigma_d))
     sigma_d = sqrt(3)*radius/4/2;
end

if(isempty(sigma_theta))
     sigma_theta = pi/(2*sqrt(3))/2;
end

if(~exist('save_tif_flag','var'))
    save_tif_flag = 0;
end

sigma_gaussian_ratio =     sigma_gaussian/((3*radius/8));
sigma_d_ratio =     sigma_d/((sqrt(3)*radius/4/2)) ;
sigma_theta_ratio =     sigma_theta/((pi/(2*sqrt(3))/2));
    

display( 'Start test for : ')
display( 'radius  sigma_gaussian   sigma_d   sigma_theta');
    
[radius  sigma_gaussian  sigma_d sigma_theta]

display( ['That is sg ',num2str(round(sigma_gaussian_ratio*100),'%d'),...
        '%,  sd',num2str(round(sigma_d_ratio*100),'%d'),...
        '%,  sa',num2str(round(sigma_theta_ratio*100),'%d'),'%']);
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
         
% outdir = [MD_1.processes_{indexFilamentSegmentationProcess_1}.outFilePaths_{iChannel1},filesep,'similarity_results_'];
% make the folder in root and with channel info

outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),'_similarity_results'];

if(~exist(outdir,'dir'))
    mkdir(outdir);
end

if(flag_default==1)
    customized_outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),...
        '_similarity_r',num2str(radius),'_s_default'];
else
%     customized_outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),...
%         '_similarity_r',num2str(radius),'_sg',num2str(sigma_gaussian_ratio,'%.2f'),...
%         '_sd',num2str(sigma_d_ratio,'%.2f'),...
%         '_sa',num2str(sigma_theta_ratio,'%.2f')];
%     customized_outdir(customized_outdir=='.')='p';
%     
    
      customized_outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),...
        '_similarity_r',num2str(radius),'_sg',num2str(round(sigma_gaussian_ratio*100),'%d'),...
        '_sd',num2str(round(sigma_d_ratio*100),'%d'),...
        '_sa',num2str(round(sigma_theta_ratio*100),'%d')];
%     customized_outdir(customized_outdir=='.')='p';
end


if(~exist(customized_outdir,'dir'))
    mkdir(customized_outdir);
end


for iFrame = start_frame1 : nFrame
    
    iFrame_2 = iFrame - start_frame1 + start_frame2;
    
    if iFrame_2<= nFrame
        disp(['Frame: ',num2str(iFrame)]);       
        
        % load first movie
        MT_orientation = MD_1.processes_{indexFilamentSegmentationProcess_1}.loadChannelOutput(iChannel1,iFrame,'output','current_seg_orientation');
        MT_current_model = MD_1.processes_{indexFilamentSegmentationProcess_1}.loadChannelOutput(iChannel1,iFrame,'output','current_model');
        
        MT_current_seg = (isnan(MT_orientation)==0);
        
        try
            MT_img =  MD_1.processes_{indexFlattenProcess_1}.loadChannelOutput(iChannel1,iFrame);
        catch
            MT_img =  MD_1.channels_(1).loadImage(iFrame);
        end
        
        %
        % load second movie, at the defined frame number
        VIF_orientation = MD_2.processes_{indexFilamentSegmentationProcess_2}.loadChannelOutput(iChannel2,iFrame_2+0,'output','current_seg_orientation');
        VIF_current_model = MD_2.processes_{indexFilamentSegmentationProcess_2}.loadChannelOutput(iChannel2,iFrame_2+0,'output','current_model');
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        
        try
            VIF_img =  MD_2.processes_{indexFlattenProcess_2}.loadChannelOutput(iChannel2,iFrame_2);
        catch
            VIF_img =  MD.channels_(2).loadImage(iFrame);
        end
        
        % get rid of out of cell filament, if not having done so previously
            
        
        try
            MT_cell_mask =  MD_1.processes_{indexCellRefineProcess_1}.loadChannelOutput(iChannel1,iFrame);            
            MT_current_model =  ROIed_digital_filament_model(MT_current_model,size(VIF_img), MT_cell_mask);
            MT_current_seg = filament_model_to_seg_bwim(MT_current_model,size(VIF_img),[]);
        end
        
        try
            VIF_cell_mask =  MD_2.processes_{indexCellRefineProcess_2}.loadChannelOutput(iChannel2,iFrame);            
            VIF_current_model =  ROIed_digital_filament_model(VIF_current_model,size(VIF_img), VIF_cell_mask);
            VIF_current_seg = filament_model_to_seg_bwim(VIF_current_model,size(VIF_img),[]);
        end
        
        % when the image is uint16, it is usually a original image, so
        % normalize it for visualizetion
        
        if( isa(VIF_img,'uint16'))
            VIF_img = double(VIF_img);          
            VIF_img = uint8(255*(VIF_img)/(max(max(VIF_img))));
        end
        
        if( isa(MT_img,'uint16'))
            MT_img = double(MT_img);
            MT_img = uint8(255*MT_img/(max(max(MT_img))));
        end
        
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
        
%         [similarity_scoremap, similarity_scoremap_1to2, ...
%             similarity_scoremap_2to1,distance_map_1_2, distance_map_2_1, angle_map_1_2, angle_map_2_1] = ...
%             network_similarity_scoremap(MT_current_model,VIF_current_model,img_size, radius,...
%             outdir,iFrame,save_everything_flag);
        
%         display('Start calculation');
        tic
        [similarity_scoremap, difference_map]= ...
            network_similarity_scoremap(MT_current_model,VIF_current_model,img_size, radius,...
            longest_radius,sigma_d, sigma_theta,sigma_gaussian);
        toc
        
%         display('Start plotting');
        % plot the detailed results if requested.
        plot_differnce_map_wrapper(difference_map,outdir,customized_outdir,iFrame,radius,show_save_everything_flag,save_tif_flag);        
        
        save([outdir,filesep,'Similarity_maps_frame',num2str(iFrame),'.mat'], 'similarity_scoremap', 'difference_map');
        save([customized_outdir,filesep,'Similarity_maps_frame',num2str(iFrame),'.mat'], 'similarity_scoremap', 'difference_map');
            
        
%         [filament_similiarity_1,filament_similiarity_2] ...
%             = filament_similarity(VIF_current_model,MT_current_model,img_size, ...
%             radius,outdir,iFrame,save_everything_flag,...
%             difference_map.distance_map_1_2, difference_map.distance_map_2_1, ...
%             difference_map.angle_map_1_2, difference_map.angle_map_2_1);
%         
        
        similarity_scoremap_cell{1, iFrame} = similarity_scoremap;
        similarity_scoremap_1to2_cell{1, iFrame} =  difference_map.similarity_scoremap_1to2;
        similarity_scoremap_2to1_cell{1, iFrame} =  difference_map.similarity_scoremap_2to1;
        distance_map_1_2_cell{1, iFrame} =  difference_map.distance_map_1_2;
        distance_map_2_1_cell{1, iFrame} =  difference_map.distance_map_2_1;
        angle_map_1_2_cell{1, iFrame} =  difference_map.angle_map_1_2;
        angle_map_2_1_cell{1, iFrame} =  difference_map.angle_map_2_1;
         
        %       close all;
    end
end

        
save([outdir,filesep,'VIFMT_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','similarity_scoremap_1to2_cell',...
    'similarity_scoremap_2to1_cell', ...
    'distance_map_1_2_cell','distance_map_2_1_cell', ...
    'angle_map_1_2_cell','angle_map_2_1_cell');
        
save([customized_outdir,filesep,'VIFMT_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','similarity_scoremap_1to2_cell',...
    'similarity_scoremap_2to1_cell', ...
    'distance_map_1_2_cell','distance_map_2_1_cell', ...
    'angle_map_1_2_cell','angle_map_2_1_cell');


% winopen(outdir);
