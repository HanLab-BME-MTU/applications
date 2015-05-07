function [similarity_scoremap_cell,difference_map_cell,similarity_scoremap_1to2_cell,similarity_scoremap_2to1_cell, distance_map_1_2_cell, distance_map_2_1_cell, angle_map_1_2_cell, angle_map_2_1_cell] ...
    = load_HFF_MD_network_for_dynamics_compare(MD1_filename,iChannel1, start_frame1,...
    iChannel2, start_frame2, ...
    radius,show_save_everything_flag,...
    longest_radius,sigma_gaussian, sigma_d, sigma_theta)
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
    longest_radius = 100;
end

flag_default = 0;

if(nargin<10)
    sigma_gaussian = 3*radius/8;
    flag_default = 1;
end
if(nargin<11)
    sigma_d = sqrt(3)*radius/4;
end

if(nargin<12)
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
    

% input: save_everything_flag:
MD = load(MD1_filename);
MD = MD.MD;

indexVIFChannel=2;
indexMTChannel=1;


package_process_ind_script;
indexFilamentSegmentationProcess = indexFilamentSegmentationProcess;
indexCellRefineProcess = indexCellRefineProcess;
indexCellSegProcess = indexCellSegProcess;


nFrame = MD.nFrames_;

%initialize output
similarity_scoremap_cell = cell(1,nFrame);
difference_map_cell = cell(1,nFrame);
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
      customized_outdir = [MD.outputDirectory_,filesep,'ch',num2str(iChannel1),'_',num2str(iChannel2),...
        '_similarity_r',num2str(radius),'_sg',num2str(round(sigma_gaussian_ratio*100),'%d'),...
        '_sd',num2str(round(sigma_d_ratio*100),'%d'),...
        '_sa',num2str(round(sigma_theta_ratio*100),'%d')];
end

if(~exist(customized_outdir,'dir'))
    mkdir(customized_outdir);
end

for iFrame = start_frame1 : nFrame
    
    iFrame_2 = iFrame - start_frame1 + start_frame2;
    
    if iFrame_2<= nFrame
        disp(['Frame: ',num2str(iFrame)]);       
        
        % load first movie
        MT_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel1,iFrame,'output','current_seg_orientation');
        MT_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel1,iFrame,'output','current_model');
        
        MT_current_seg = (isnan(MT_orientation)==0);

        %
        % load second movie, at the defined frame number
        VIF_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel2,iFrame_2+0,'output','current_seg_orientation');
        VIF_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel2,iFrame_2+0,'output','current_model');
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        
        img_size = size(MT_orientation);
        
        tic
        [similarity_scoremap, difference_map]= ...
            network_similarity_scoremap(MT_current_model,VIF_current_model,img_size, radius,...
            longest_radius,sigma_d, sigma_theta,sigma_gaussian);
        toc
        
%         display('Start plotting');
        % plot the detailed results if requested.
        plot_differnce_map_wrapper(difference_map,outdir,customized_outdir,iFrame,radius,show_save_everything_flag);        
        
        save([outdir,filesep,'Similarity_maps_frame',num2str(iFrame),'.mat'], 'similarity_scoremap', 'difference_map');
        save([customized_outdir,filesep,'Similarity_maps_frame',num2str(iFrame),'.mat'], 'similarity_scoremap', 'difference_map');
            
        difference_map_cell{1, iFrame} = difference_map;
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

        
save([outdir,filesep,'VIF_Actin_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','similarity_scoremap_1to2_cell',...
    'similarity_scoremap_2to1_cell', ...
    'distance_map_1_2_cell','distance_map_2_1_cell', ...
    'angle_map_1_2_cell','angle_map_2_1_cell','difference_map_cell');
        
save([customized_outdir,filesep,'VIFMT_Actin_sm_maps_allframe.mat'], ...
    'similarity_scoremap_cell','similarity_scoremap_1to2_cell',...
    'similarity_scoremap_2to1_cell', ...
    'distance_map_1_2_cell','distance_map_2_1_cell', ...
    'angle_map_1_2_cell','angle_map_2_1_cell','difference_map_cell');


% winopen(outdir);
