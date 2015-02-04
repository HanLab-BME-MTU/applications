function NA_Group = network_analysis_group_plotting(ML_name_cell, feature_index, Group_ROOT_DIR)
% function to do network feature plotting for a whole movielist and then
% for a list of ML
% Liya Ding, March, 2014
%
% Input:
%   ML_array:     The array of movieList objects in this group to be analyzed


if(nargin<2)
    feature_index=ones(18,1);
end

if(nargin<3)
    Group_ROOT_DIR=[];
end

close all;

nList = length(ML_name_cell);

for iML = 1 : nList
    
   try
       ML = MovieList.load(ML_name_cell{iML});
   catch
       display('Error in ML loading.');
       break;
   end
   
   ML_ROOT_DIR = ML.outputDirectory_;
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    NA_feature_thisML = cell(1,1);
    Identifier_thisML = cell(1,1);
    
    % if the batch results exist, load it
    if 0
        
        % this part commented because the result for a whole ML will be too
        % large to store
%         (exist([ML_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'], 'file'))
%         load([ML_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'],...
%             'NA_feature_thisML','Identifier_thisML');
    else
        for iMD  = 1 : movieNumber
            MD_ROOT_DIR = MD.outputDirectory_;
            
            % otherwise load one by one
            NA_feature_thisMD = cell(1,1);
            Identifier_thisMD = cell(1,1);
            
            if(exist([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'], 'file'))
                load([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'],...
                    'NA_feature_thisMD','Identifier_thisMD');
            else                
                try
                    % load this movie
                    MD =  MovieData.load(ML.movieDataFile_{iMD});
                catch
                    disp('The MD file is missing');
                    continue;
                end
                
                % Now MD in workspace
                %check index
                package_process_ind_script;
                
                nChannel = numel(MD.channels_);
                nFrame = MD.nFrames_;
                
                for iChannel = 1 : nChannel
                    display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);
                    outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
                    
                    Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
                    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
                    
                    % check if this folder exist
                    if(exist(outdir,'dir'))
                        % if it exist, try to do the branch analysis
                        %                         try
                        
                        for iFrame = 1 : nFrame
                            display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
                            load_ch_frame_flag =0;
                            
                            if(exist([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),'_',...
                                    '_frame_',num2str(iFrame),'.mat'],'file'))
                                load([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),...
                                    '_frame_',num2str(iFrame),'.mat'],'output_feature');
                                load_ch_frame_flag=1;
                            end
                            
                            if(exist([outdir,filesep,'network_feature_ch',num2str(iChannel),'_',...
                                    filename_short_strs{iFrame},'.mat'],'file'))
                                load([outdir,filesep,'network_feature_ch',num2str(iChannel),'_',...
                                    filename_short_strs{iFrame},'.mat'], 'output_feature');
                                load_ch_frame_flag = 1;
                            end
                            
                            if(load_ch_frame_flag==0)
                                continue;
                            end
                            
                            display(['Found: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
                            Identifier = ['ML',num2str(iML),'M',num2str(iMD),'F',num2str(iFrame)];
                            
                            % save
                            NA_feature_thisMD{iChannel, iFrame} = output_feature;
                            Identifier_thisMD{iChannel, iFrame} = filename_short_strs{iFrame};                            
                            CFMP_feature_thisMD{iChannel, iFrame} = extract_mean_percentiles_one(output_feature, requested_feature_index);
                            
                            %                         end
                            
                        end% end of a frame
                        
                    end % end of if the analysis for output folder exist
                    
                    ChMP_feature_thisMD{1, iFrame} = extract_mean_percentiles_cell(NA_feature_thisMD{1, iFrame}, requested_feature_index);
                                  
                    
                end % end of a Channel
                
                % for all informaion in this MD, gather into cell
                NA_feature_MD_cell{iMD} = NA_feature_thisMD;
                Identifier_MD_cell{iMD}= Identifier_thisMD;
                
                CFMP_feature_MD_cell{iMD}= CFMP_feature_thisMD;
                ChMP_feature_MD_cell{iMD}= ChMP_feature_thisMD;
                
                % save for results this whole movie
                save([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'],...
                    'NA_feature_thisMD','Identifier_thisMD','CFMP_feature_MD_cell','ChMP_feature_MD_cell');
                
            end % end of if previous gathering exists
            
            NA_feature_thisML{1,iMD} = NA_feature_thisMD;
            Identifier_thisML{1,iMD} = Identifier_thisMD;
            CFMP_feature_thisML{1,iMD} = CFMP_feature_thisMD;
            ChMP_feature_thisML{1,iMD} = ChMP_feature_thisMD;
        
        end  % end of a MD        
        
       
    end % end of a ML    
    
    
end

% initialize the pools

Group_Pool_Travel_Length = [];
Group_Pool_Travel_Distance = [];
Group_Pool_Travel_Speed = [];
Group_Pool_Cell_Marked_Frame_Number = [];
Group_Pool_branch_number_tracked = [];
Group_Pool_branch_number_max = [];
Group_Pool_branch_number_mean = [];
Group_Pool_branch_duration_array = [];
Group_Pool_branch_duration_mean= [];
Group_Pool_branch_vif_mean_intensity= [];
Group_Pool_protrusion_vif_mean_intensity= [];
Group_Pool_retraction_vif_mean_intensity= [];
Group_Pool_whole_cell_vif_mean_intensity= [];
Group_Pool_whole_cell_size_mean= [];
Group_Pool_branch_number_mean_pat=[];
Group_Pool_whole_cell_vim_totalamount_mean=[];
Group_Pool_whole_cell_vif_mean_intensity_pat=[];
Group_Pool_branch_size_mean=[];
Group_Pool_thresholded_branch_number_mean=[];
Group_Pool_branch_cellmovement_std=[];
Group_Pool_fila_branch_orientation_pool_std=[];
Group_Pool_fila_trajectory_orientation_pool_std=[];
Group_Pool_fila_branch_orientation_pool=[];
Group_Pool_fila_trajectory_orientation_pool=[];
Group_Pool_fila_branch_orientation_pool_slow=[];
Group_Pool_fila_trajectory_orientation_pool_slow=[];
Group_Pool_fila_branch_orientation_pool_fast=[];
Group_Pool_fila_trajectory_orientation_pool_fast=[];
Group_Pool_pool_vif_intensity=[];
Group_Pool_CompletedFrame_last=[];
Group_Pool_whole_cell_vim_seg_total = [];
Group_Pool_whole_cell_vim_seg_mean  = [];
Group_Pool_whole_cell_vim_nms_total = [];
Group_Pool_whole_cell_vim_nms_mean  = [];
Group_Pool_branch_seg_total = [];
Group_Pool_branch_seg_mean  = [];
Group_Pool_branch_nms_total = [];
Group_Pool_branch_nms_mean  = [];


Identifier_cell = [];

iAllCell = 0;

% if no input, the DIR of last movieList
if(isempty(Group_ROOT_DIR))
    Group_ROOT_DIR = ML_ROOT_DIR;
end

if(exist(Group_ROOT_DIR,'dir')==0)
    mkdir(Group_ROOT_DIR);
end

%% For each ML
for iML = 1 : nList
    NA_feature_ML_thisMD = NA_Group{1, iML};
    
    ML_Pool_Travel_Length = [];
    ML_Pool_Travel_Distance = [];
    ML_Pool_Travel_Speed = [];
    ML_Pool_Cell_Marked_Frame_Number = [];
    ML_Pool_branch_number_tracked = [];
    ML_Pool_branch_number_max = [];
    ML_Pool_branch_number_mean = [];
    ML_Pool_branch_duration_array = [];
    ML_Pool_branch_duration_mean= [];
    ML_Pool_branch_vif_mean_intensity= [];
    ML_Pool_protrusion_vif_mean_intensity= [];
    ML_Pool_retraction_vif_mean_intensity= [];
    ML_Pool_whole_cell_vif_mean_intensity= [];
    ML_Pool_whole_cell_size_mean= [];
    ML_Pool_branch_number_mean_pat=[];
    ML_Pool_whole_cell_vim_totalamount_mean=[];
    ML_Pool_whole_cell_vif_mean_intensity_pat=[];
    ML_Pool_branch_size_mean=[];
    ML_Pool_thresholded_branch_number_mean=[];
    ML_Pool_branch_cellmovement_std=[];
    ML_Pool_fila_branch_orientation_pool_std=[];
    ML_Pool_fila_trajectory_orientation_pool_std=[];
    ML_Pool_fila_branch_orientation_pool=[];
    ML_Pool_fila_trajectory_orientation_pool=[];
    ML_Pool_fila_branch_orientation_pool_slow=[];
    ML_Pool_fila_trajectory_orientation_pool_slow=[];
    ML_Pool_fila_branch_orientation_pool_fast=[];
    ML_Pool_fila_trajectory_orientation_pool_fast=[];
    ML_Pool_CompletedFrame_last = [];
    ML_Pool_pool_vif_intensity=[];
    ML_Pool_whole_cell_vim_seg_total = [];
    ML_Pool_whole_cell_vim_seg_mean = [];
    ML_Pool_whole_cell_vim_nms_total = [];
    ML_Pool_whole_cell_vim_nms_mean = [];
    ML_Pool_branch_seg_total = [];
    ML_Pool_branch_seg_mean = [];
    ML_Pool_branch_nms_total = [];
    ML_Pool_branch_nms_mean = [];
    
    % the number of movies
    movieNumber =  20;
    
    for iMD  = 1 : movieNumber
        nChannel=2;
        for iChannel = 1 : nChannel
            for iFrame = 1 :  10
                try
                    NA_feature = NA_feature_ML_thisMD{1, iMD}{iChannel, iFrame};
                catch
                    continue;
                end
                
                if(~isempty(NA_feature))
                    ML_Pool_CompletedFrame_last=[ML_Pool_CompletedFrame_last NA_feature.CompletedFrame(end)];
                    
                    ML_Pool_branch_size_mean = [ML_Pool_branch_size_mean NA_feature.branch_mean_size];
                    ML_Pool_whole_cell_size_mean= [ML_Pool_whole_cell_size_mean NA_feature.whole_cell_size_mean];
                    ML_Pool_branch_number_mean_pat = [ML_Pool_branch_number_mean_pat repmat(NA_feature.branch_number_mean, [1, length(NA_feature.branch_vif_mean_intensity)])];
                    Identifier_cell{1,numel(Identifier_cell)+1} = NA_feature.Identifier;
                    
                    
                    Thresholded_branch_number_mean = sum(NA_feature.branch_duration_array(find(NA_feature.branch_mean_size>T_branchsize & NA_feature.branch_duration_array>T_branchduration)))...
                        ./ NA_feature.cell_marked_frame_number;
                    ML_Pool_thresholded_branch_number_mean = [ML_Pool_thresholded_branch_number_mean Thresholded_branch_number_mean];
                    
                    
                    ML_Pool_Travel_Length = [ML_Pool_Travel_Length NA_feature.cell_travel_length];
                    ML_Pool_Travel_Distance = [ML_Pool_Travel_Distance NA_feature.cell_travel_distance];
                    ML_Pool_Travel_Speed = [ML_Pool_Travel_Speed NA_feature.cell_travel_speed];
                    ML_Pool_Cell_Marked_Frame_Number = [ML_Pool_Cell_Marked_Frame_Number NA_feature.cell_marked_frame_number];
                    ML_Pool_branch_number_tracked = [ML_Pool_branch_number_tracked NA_feature.branch_number_tracked];
                    ML_Pool_branch_number_max = [ML_Pool_branch_number_max NA_feature.branch_number_max];
                    ML_Pool_branch_number_mean = [ML_Pool_branch_number_mean NA_feature.branch_number_mean];
                    ML_Pool_branch_duration_array = [ML_Pool_branch_duration_array  NA_feature.branch_duration_array];
                    ML_Pool_branch_duration_mean= [ML_Pool_branch_duration_mean  NA_feature.branch_duration_mean];
                    
                    ML_Pool_whole_cell_vif_mean_intensity_pat = [ML_Pool_whole_cell_vif_mean_intensity_pat repmat(NA_feature.whole_cell_vif_mean_intensity, [1, length(NA_feature.branch_vif_mean_intensity)])];
                    ML_Pool_whole_cell_vim_totalamount_mean = [ML_Pool_whole_cell_vim_totalamount_mean NA_feature.whole_cell_vim_totalamount_mean];
                    ML_Pool_branch_vif_mean_intensity= [ML_Pool_branch_vif_mean_intensity NA_feature.branch_vif_mean_intensity];
                    ML_Pool_protrusion_vif_mean_intensity= [ML_Pool_protrusion_vif_mean_intensity NA_feature.protrusion_vif_mean_intensity];
                    ML_Pool_retraction_vif_mean_intensity= [ML_Pool_retraction_vif_mean_intensity NA_feature.retraction_vif_mean_intensity];
                    ML_Pool_whole_cell_vif_mean_intensity= [ML_Pool_whole_cell_vif_mean_intensity NA_feature.whole_cell_vif_mean_intensity];
                    
                    ML_Pool_pool_vif_intensity= [ML_Pool_pool_vif_intensity; NA_feature.pool_all_vif_intensity];
                    
                    
                    
                    ML_Pool_whole_cell_vim_seg_total= [ML_Pool_whole_cell_vim_seg_total; NA_feature.whole_cell_vim_seg_total];
                    ML_Pool_whole_cell_vim_seg_mean= [ML_Pool_whole_cell_vim_seg_mean; NA_feature.whole_cell_vim_seg_mean];
                    ML_Pool_whole_cell_vim_nms_total= [ML_Pool_whole_cell_vim_nms_total; NA_feature.whole_cell_vim_nms_total];
                    ML_Pool_whole_cell_vim_nms_mean= [ML_Pool_whole_cell_vim_nms_mean; NA_feature.whole_cell_vim_nms_mean];
                    
                    ML_Pool_branch_seg_total= [ML_Pool_branch_seg_total NA_feature.branch_seg_total];
                    ML_Pool_branch_seg_mean= [ML_Pool_branch_seg_mean NA_feature.branch_seg_mean];
                    ML_Pool_branch_nms_total= [ML_Pool_branch_nms_total NA_feature.branch_nms_total];
                    ML_Pool_branch_nms_mean= [ML_Pool_branch_nms_mean NA_feature.branch_nms_mean];
                    
                    
                    
                    % some old analysis doens't have orient std,
                    % but as space holder, add nan
                    try
                        ML_Pool_branch_cellmovement_std= [ML_Pool_branch_cellmovement_std NA_feature.branch_cellmovement_std];
                        
                    catch
                        ML_Pool_branch_cellmovement_std= [ML_Pool_branch_cellmovement_std nan];
                    end
                    
                    % some old analysis doens't have orient std,
                    % but as space holder, add nan
                    try
                        ML_Pool_fila_branch_orientation_pool_std= [ML_Pool_fila_branch_orientation_pool_std NA_feature.fila_branch_orientation_pool_std];
                        
                    catch
                        ML_Pool_fila_branch_orientation_pool_std= [ML_Pool_fila_branch_orientation_pool_std nan];
                    end
                    
                    try
                        ML_Pool_fila_trajectory_orientation_pool_std= [ML_Pool_fila_trajectory_orientation_pool_std NA_feature.fila_trajectory_orientation_pool_std];
                        
                    catch
                        ML_Pool_fila_trajectory_orientation_pool_std= [ML_Pool_fila_trajectory_orientation_pool_std nan];
                    end
                    
                    %some old analysis doesn't have these pools
                    try
                        ML_Pool_branch_cellmovement_pool= [ML_Pool_branch_cellmovement_pool; NA_feature.branch_trajectory_orientation_pool];
                        
                        if(NA_feature.cell_travel_speed>7)
                            ML_Pool_branch_cellmovement_pool_fast= [ML_Pool_branch_cellmovement_pool_fast; NA_feature.branch_trajectory_orientation_pool];
                        else
                            ML_Pool_branch_cellmovement_pool_slow= [ML_Pool_branch_cellmovement_pool_slow; NA_feature.branch_trajectory_orientation_pool];
                        end
                    end
                    
                    
                    
                    try
                        ML_Pool_fila_branch_orientation_pool= [ML_Pool_fila_branch_orientation_pool; NA_feature.fila_branch_orientation_pool];
                        
                        if(NA_feature.cell_travel_speed>7)
                            ML_Pool_fila_branch_orientation_pool_fast=[ML_Pool_fila_branch_orientation_pool_fast; NA_feature.fila_branch_orientation_pool];
                        else
                            ML_Pool_fila_branch_orientation_pool_slow=[ML_Pool_fila_trajectory_orientation_pool_slow; NA_feature.fila_trajectory_orientation_pool];
                        end
                    end
                    
                    
                    try
                        ML_Pool_fila_trajectory_orientation_pool= [ML_Pool_fila_trajectory_orientation_pool; NA_feature.fila_trajectory_orientation_pool];
                        
                        if(NA_feature.cell_travel_speed>7)
                            ML_Pool_fila_trajectory_orientation_pool_fast=[ML_Pool_fila_trajectory_orientation_pool_fast; NA_feature.fila_trajectory_orientation_pool];
                        else
                            ML_Pool_fila_trajectory_orientation_pool_slow=[ML_Pool_fila_branch_orientation_pool_slow; NA_feature.fila_branch_orientation_pool];
                        end
                    end
                    iAllCell = iAllCell + 1;
                end
            end
        end
    end
    
    %% vim normalization
    vim_max = prctile(double(ML_Pool_pool_vif_intensity),98);
    
    display('50')
    prctile(double(ML_Pool_pool_vif_intensity),50)
    display('75')
    prctile(double(ML_Pool_pool_vif_intensity),75)
    display('95')
    prctile(double(ML_Pool_pool_vif_intensity),95)
    display('98')
    prctile(double(ML_Pool_pool_vif_intensity),98)
    display('99')
    prctile(double(ML_Pool_pool_vif_intensity),99)
    
    
    ML_Pool_branch_vif_mean_intensity = ML_Pool_branch_vif_mean_intensity/vim_max;
    ML_Pool_protrusion_vif_mean_intensity = ML_Pool_protrusion_vif_mean_intensity/vim_max;
    ML_Pool_retraction_vif_mean_intensity = ML_Pool_retraction_vif_mean_intensity/vim_max;
    ML_Pool_whole_cell_vif_mean_intensity = ML_Pool_whole_cell_vif_mean_intensity/vim_max;
    ML_Pool_whole_cell_vif_mean_intensity_pat = ML_Pool_whole_cell_vif_mean_intensity_pat/vim_max;
    ML_Pool_whole_cell_vim_totalamount_mean = ML_Pool_whole_cell_vim_totalamount_mean/vim_max;
    
    %% in the end put data from this ML to the pool of all MLs
    
    Group_Pool_Travel_Length = [Group_Pool_Travel_Length ML_Pool_Travel_Length];
    Group_Pool_Travel_Distance = [Group_Pool_Travel_Distance ML_Pool_Travel_Distance];
    Group_Pool_Travel_Speed = [Group_Pool_Travel_Speed ML_Pool_Travel_Speed];
    Group_Pool_Cell_Marked_Frame_Number = [Group_Pool_Cell_Marked_Frame_Number ML_Pool_Cell_Marked_Frame_Number];
    Group_Pool_branch_number_tracked = [Group_Pool_branch_number_tracked ML_Pool_branch_number_tracked];
    Group_Pool_branch_number_max = [Group_Pool_branch_number_max ML_Pool_branch_number_max];
    Group_Pool_branch_number_mean = [Group_Pool_branch_number_mean ML_Pool_branch_number_mean];
    Group_Pool_branch_duration_array = [Group_Pool_branch_duration_array ML_Pool_branch_duration_array];
    Group_Pool_branch_duration_mean= [Group_Pool_branch_duration_mean ML_Pool_branch_duration_mean];
    Group_Pool_branch_vif_mean_intensity= [Group_Pool_branch_vif_mean_intensity ML_Pool_branch_vif_mean_intensity];
    Group_Pool_protrusion_vif_mean_intensity= [Group_Pool_protrusion_vif_mean_intensity ML_Pool_protrusion_vif_mean_intensity];
    Group_Pool_retraction_vif_mean_intensity= [Group_Pool_retraction_vif_mean_intensity ML_Pool_retraction_vif_mean_intensity];
    Group_Pool_whole_cell_vif_mean_intensity= [Group_Pool_whole_cell_vif_mean_intensity ML_Pool_whole_cell_vif_mean_intensity];
    Group_Pool_whole_cell_size_mean= [Group_Pool_whole_cell_size_mean ML_Pool_whole_cell_size_mean];
    Group_Pool_branch_number_mean_pat=[Group_Pool_branch_number_mean_pat ML_Pool_branch_number_mean_pat];
    Group_Pool_whole_cell_vim_totalamount_mean=[Group_Pool_whole_cell_vim_totalamount_mean ML_Pool_whole_cell_vim_totalamount_mean];
    Group_Pool_whole_cell_vif_mean_intensity_pat=[Group_Pool_whole_cell_vif_mean_intensity_pat ML_Pool_whole_cell_vif_mean_intensity_pat];
    Group_Pool_branch_size_mean=[Group_Pool_branch_size_mean ML_Pool_branch_size_mean];
    Group_Pool_thresholded_branch_number_mean=[Group_Pool_thresholded_branch_number_mean ML_Pool_thresholded_branch_number_mean];
    Group_Pool_branch_cellmovement_std=[Group_Pool_branch_cellmovement_std ML_Pool_branch_cellmovement_std];
    Group_Pool_fila_branch_orientation_pool_std=[Group_Pool_fila_branch_orientation_pool_std ML_Pool_fila_branch_orientation_pool_std];
    Group_Pool_fila_trajectory_orientation_pool_std=[Group_Pool_fila_trajectory_orientation_pool_std ML_Pool_fila_trajectory_orientation_pool_std];
    Group_Pool_fila_branch_orientation_pool=[Group_Pool_fila_branch_orientation_pool; ML_Pool_fila_branch_orientation_pool];
    Group_Pool_fila_trajectory_orientation_pool=[Group_Pool_fila_trajectory_orientation_pool; ML_Pool_fila_trajectory_orientation_pool];
    Group_Pool_fila_branch_orientation_pool_slow=[Group_Pool_fila_branch_orientation_pool_slow; ML_Pool_fila_branch_orientation_pool_slow];
    Group_Pool_fila_trajectory_orientation_pool_slow=[Group_Pool_fila_trajectory_orientation_pool_slow; ML_Pool_fila_trajectory_orientation_pool_slow];
    Group_Pool_fila_branch_orientation_pool_fast=[Group_Pool_fila_branch_orientation_pool_fast; ML_Pool_fila_branch_orientation_pool_fast];
    Group_Pool_fila_trajectory_orientation_pool_fast=[Group_Pool_fila_trajectory_orientation_pool_fast; ML_Pool_fila_trajectory_orientation_pool_fast];
    Group_Pool_CompletedFrame_last =[Group_Pool_CompletedFrame_last ML_Pool_CompletedFrame_last];
    
    
    Group_Pool_whole_cell_vim_seg_total = [Group_Pool_whole_cell_vim_seg_total; ML_Pool_whole_cell_vim_seg_total];
    Group_Pool_whole_cell_vim_seg_mean  = [Group_Pool_whole_cell_vim_seg_mean;  ML_Pool_whole_cell_vim_seg_mean ];
    Group_Pool_whole_cell_vim_nms_total = [Group_Pool_whole_cell_vim_nms_total; ML_Pool_whole_cell_vim_nms_total];
    Group_Pool_whole_cell_vim_nms_mean  = [Group_Pool_whole_cell_vim_nms_mean; ML_Pool_whole_cell_vim_nms_mean ];
    
    Group_Pool_branch_seg_total = [Group_Pool_branch_seg_total ML_Pool_branch_seg_total];
    Group_Pool_branch_seg_mean  = [Group_Pool_branch_seg_mean  ML_Pool_branch_seg_mean ];
    Group_Pool_branch_nms_total = [Group_Pool_branch_nms_total ML_Pool_branch_nms_total];
    Group_Pool_branch_nms_mean  = [Group_Pool_branch_nms_mean  ML_Pool_branch_nms_mean ];
    
end

%% Plotting for all MLs
Group_Pool_tracked_branch_d_frames = Group_Pool_branch_number_tracked./Group_Pool_Cell_Marked_Frame_Number;

colorarray = rand(numel(Group_Pool_branch_vif_mean_intensity),3);



%% % plot the this group result
h3 = figure(3);
plot(Group_Pool_branch_duration_array, Group_Pool_branch_vif_mean_intensity,'.');
xlabel(['Branch Duration,sample size:',num2str(numel(Group_Pool_branch_duration_array))],'Fontsize',13);
ylabel('Branch Vim mean Int','Fontsize',13);
title('Branch: Duration vs Vim Mean Int','Fontsize',13);
set(gca,'fontsize',13)
saveas(h3,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.fig']);
saveas(h3,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.tif']);
print(h3,'-depsc',[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.eps']);


size(Group_Pool_protrusion_vif_mean_intensity)

h4 = figure(4);
[h_p,bin] = hist(Group_Pool_protrusion_vif_mean_intensity,0.0:0.03:1.1);
[h_r,bin] = hist(Group_Pool_retraction_vif_mean_intensity,0.0:0.03:1.1);
bar([reshape(bin, 1, numel(bin)); reshape(bin, 1, numel(bin))]',[reshape(h_p, 1, numel(h_p));reshape(h_r, 1, numel(h_r))]');
legend('Protrusion Region Vim','Retraction Region Vim');
axis([0.00 0.30 0 max([h_p h_r])])
set(gca,'fontsize',13)
saveas(h4,[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.fig']);
saveas(h4,[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.tif']);
print(h4,'-depsc',[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.eps']);

% mean vim vs speed
h7 = figure(7); hold off;
%  plot(Group_Pool_branch_number_mean, Group_Pool_Travel_Speed,'o');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_Travel_Speed(i)+0.2, num2str(i), 'color',colorarray(i,:))
    % hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'o', 'color',colorarray(i,:))
    hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'bo','linewidth',2,'markersize',7);
    
end
xlabel(['Average Branch Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'FontSize',13);
ylabel('Cell Speed','FontSize',13);
title(['Branchness vs Speed, ',...
    ' Correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean', ...
    Group_Pool_Travel_Speed'),'%1.2f')],'FontSize',13);
% axis([0 10 0 24]);
set(gca,'fontsize',13);
saveas(h7,[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.fig']);
saveas(h7,[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.tif']);
print(h7,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.eps']);


h8 = figure(8);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16)
    plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
end
xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'Fontsize',13);
ylabel('Cell Vim mean Int','Fontsize',13);
title(['Branch number vs Vim, correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
set(gca,'fontsize',13);
saveas(h8,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.fig']);
saveas(h8,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.tif']);
print(h8,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.eps']);


h8_axis = axis;

h18 = figure(18);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    %     text(Group_Pool_thresholded_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
end
corr_b_v_v = corrcoef(Group_Pool_whole_cell_vif_mean_intensity',Group_Pool_thresholded_branch_number_mean');
display(['Corr for Branchness vs Vim: ',num2str(corr_b_v_v(1,2))]);

xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_thresholded_branch_number_mean))],'Fontsize',13);
ylabel('Cell Vim mean Int','Fontsize',13);
title({['Branch number vs Vim, for branches larger than ',num2str(T_branchsize),' pixels'],['duration longer than ',num2str(T_branchduration),' frames'],...
    ['Correlation:',num2str(corr_b_v_v(1,2))]},'Fontsize',13);
axis(h8_axis);
set(gca,'fontsize',13);
saveas(h18,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.fig']);
saveas(h18,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.tif']);
print(h18,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.eps']);



%%
try
    %%
    h26 = figure(26);hold off;
    % plot(Group_Pool_fila_branch_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
        % text(Group_Pool_fila_branch_orientation_pool_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
        hold on;
        % plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
        
    end
    % axis([0.55, 1.0, 350 950]);
    xlabel('Filament Orientation along Branch, Standard Deviation','Fontsize',13);
    ylabel('Cell Vim mean Int','Fontsize',13);
    title({'Filament Orientation Scatterness vs Vim Mean Int',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
        Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f'), ...
        ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
    set(gca,'fontsize',13);
    saveas(h26,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.fig']);
    saveas(h26,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.tif']);
    
    print(h26,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.eps']);
    
    
    %%
    
    h36 = figure(36);hold off;
    % plot(Group_Pool_fila_branch_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_branch_number_mean)
        % text(Group_Pool_fila_branch_orientation_pool_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
        hold on;
        %         plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_branch_number_mean(i), 'bo','linewidth',2,'markersize',7);
    end
    xlabel('Filament Orientation along Branch, Standard Deviation','Fontsize',13);
    ylabel('Branch Mean Number ','Fontsize',13);
    title({'Filament Orient along Branch Scatterness vs Branch Mean Num',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))', ...
        Group_Pool_branch_number_mean(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f'),...
        ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))))]},'Fontsize',13);
    axis([0.6 0.9 0 10]);
    set(gca,'fontsize',13);
    
    saveas(h36,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.fig']);
    saveas(h36,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.tif']);
    print(h36,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.eps']);
    
    
end

try
    h46 = figure(46);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    
    title('Filament Orientation along Branch Orientaion Distribution','Fontsize',13,'fontname','Helvetica');
    set(gca,'fontsize',13);
    saveas(h46,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.fig']);
    saveas(h46,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.tif']);
    
    print(h46,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.eps']);
    
    
    
    h56 = figure(56);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool_fast,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    title('Filament Orientation along Branch Orientaion Distribution(fast cell)','Fontsize',13,'fontname','Helvetica');
    set(gca,'fontsize',13);
    saveas(h56,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.fig']);
    saveas(h56,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.tif']);
    print(h56,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.eps']);
    
    
    h66 = figure(66);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool_slow,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    set(gca,'fontsize',13);
    title('Filament Orientation along Branch Orientaion Distribution(slow cell)','Fontsize',13,'fontname','Helvetica');
    
    saveas(h66,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.fig']);
    saveas(h66,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.tif']);
    print(h66,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.eps']);
    
    
    %%
    
    %%
    h27 = figure(27);hold off;
    % plot(Group_Pool_fila_trajectory_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
        % text(Group_Pool_fila_trajectory_orientation_pool_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
        hold on;
        % plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
        
    end
    % axis([0.55, 1.0, 350 950]);
    xlabel('Filament Orientation along Path, Standard Deviation','Fontsize',13);
    ylabel('Cell Vim mean Int','Fontsize',13);
    title({'Filament Orientation along Path Scatterness vs Vim Mean Int',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
        Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f'),...
        ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
    set(gca,'fontsize',13);
    saveas(h27,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.fig']);
    saveas(h27,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.tif']);
    print(h27,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.eps']);
    
    
    %%
end

try
    h37 = figure(37);hold off;
    % plot(Group_Pool_fila_trajectory_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
    for i = 1 : length(Group_Pool_branch_number_mean)
        % text(Group_Pool_fila_trajectory_orientation_pool_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
        hold on;
        %         plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_branch_number_mean(i),...
            'bo','linewidth',2,'markersize',7);
    end
    xlabel('Filament Orientation along Cell Movement, Standard Deviation','Fontsize',13);
    ylabel('Branch Mean Number ','Fontsize',13);
    title({'Filament Orient along Cell Path Scatterness vs Branch Mean Num',...
        ['Correlation: ',...
        num2str(corr(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))', ...
        Group_Pool_branch_number_mean(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f') ,...
        ', Sample Size:',num2str(numel(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))))]},'Fontsize',13);
    set(gca,'fontsize',13);
    saveas(h37,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.fig']);
    saveas(h37,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.tif']);
    
    print(h37,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.eps']);
    
    
end

try
    h57 = figure(57);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    
    title('Filament Orientation along Cell Movement Distribution','Fontsize',13,'fontname','Helvetica');
    
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.fig']);
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.tif']);
    
    print(h57,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.eps']);
    
    h57 = figure(57);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool_fast,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    title('Filament Orientation along Cell Movement Distribution(fast cell)','Fontsize',13,'fontname','Helvetica');
    
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.fig']);
    saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.tif']);
    
    print(h57,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.eps']);
    
    
    h67 = figure(67);hold off;
    [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool_slow,-pi/2+pi/36:pi/18:pi/2+pi/36);
    bar(h_result/sum(h_result));
    
    set(gca,...
        'xlim',[0.5 19.5],...
        'xtick',[0.5:9:19.5],...
        'ylim',[0.0 0.15],...
        'ytick',[0.0:0.05:0.15],...
        'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
        'xticklabel',{'-p/2'  '0' 'p/2' },...
        'fontname','symbol',...
        'fontsize',13)
    title('Filament Orientation along Cell Movement Distribution(slow cells)','Fontsize',13,'fontname','Helvetica');
    set(gca,'fontsize',13);
    saveas(h67,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.fig']);
    saveas(h67,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.tif']);
    print(h67,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.eps']);
    
end
%%
h28 = figure(28);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    % plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
    
end
% axis([0.55, 1.0, 350 950]);
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
ylabel('Cell Vim mean Int','Fontsize',13);
title({'Branch Orientation Scatterness vs Vim',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
    Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f') ,...
    ', Sample Size:',num2str(numel(Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_whole_cell_vif_mean_intensity)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h28,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.fig']);
saveas(h28,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.tif']);
print(h28,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.eps']);


%%

h38 = figure(38);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_branch_number_mean)
    % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
ylabel('Branch Mean Number ','Fontsize',13);
title({'Branch Orientation Scatterness vs Branch Mean Number',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_branch_number_mean))', ...
    Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f') ,...
    ', Sample Size:',num2str(numel(Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_branch_cellmovement_std))))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h38,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.fig']);
saveas(h38,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.tif']);
print(h38,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.eps']);

h48 = figure(48);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : numel(Group_Pool_whole_cell_vif_mean_intensity)
    text(Group_Pool_Travel_Speed(i)-0.02, ...
        Group_Pool_whole_cell_vif_mean_intensity(i)+0.02, ...
        Identifier_cell{i}, 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Cell Speed','Fontsize',13);
ylabel('Vim Level ','Fontsize',13);
title({'Cell Speed vs Vim Level',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_Travel_Speed', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f'),...
    ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h48,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.fig']);
saveas(h48,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.tif']);
print(h48,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.eps']);



h58 = figure(58);hold off;

% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_Travel_Speed(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+0.2, num2str(i), 'color',colorarray(i,:));
    if(Group_Pool_Cell_Marked_Frame_Number(i)>15 && Group_Pool_CompletedFrame_last(i)>60)
        hold on;
        %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
        plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
    end
end
xlabel('Cell Speed','Fontsize',13);
ylabel('Vim Level ','Fontsize',13);
title({'Cell Speed vs Vim Level (cell with more frames, with late hours)',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)', ...
    Group_Pool_whole_cell_vif_mean_intensity(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)'),'%1.2f'),...
    ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h58,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.fig']);
saveas(h58,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.tif']);
print(h58,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.eps']);



%%
%
h9 = figure(9);
plot(Group_Pool_branch_number_mean_pat, Group_Pool_branch_vif_mean_intensity,'.');
xlabel('Branch Mean Number','fontsize',13);
ylabel('Branch Vim mean Int','fontsize',13);
title('Branch number vs Vim','fontsize',13);

set(gca,'fontsize',13);
saveas(h9,[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.fig']);
saveas(h9,[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.tif']);
print(h9,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.eps']);


h10 = figure(10);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vim_totalamount_mean)
    %     text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
    plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), 'bo','linewidth',2,'markersize',7);
    
end
xlabel('Branch Mean Number','fontsize',13);
ylabel('Cell Vim Total Int','fontsize',13);
title({['Branch number vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')'],...
    ['Correlation: ',num2str(corr(Group_Pool_branch_number_mean',Group_Pool_whole_cell_vim_totalamount_mean'))]},'fontsize',13);
set(gca,'fontsize',13);

saveas(h10,[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.fig']);
saveas(h10,[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.tif']);
print(h10,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.eps']);


% h20 = figure(20);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% hold on;
% plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% end
% xlabel('Branch Mean Number','Fontsize',13);
% ylabel('Cell Vim Mean Int','Fontsize',13);
% title(['Branch number vs Vim Mean,  correlation: ',...
%     num2str(corr(Group_Pool_branch_number_mean', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% saveas(h20,[Group_ROOT_DIR,filesep,'Branchness_vs_VimMean.fig']);
% saveas(h20,[Group_ROOT_DIR,filesep,'Branchness_vs_VimMean.tif']);
%
% h30 = figure(30);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% hold on;
% plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% end
% xlabel('Cell Speed','Fontsize',13);
% ylabel('Cell Vim Mean Int','Fontsize',13);
% title(['Cell Speed vs Vim Mean,  correlation: ',...
%     num2str(corr(Group_Pool_Travel_Speed', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% saveas(h30,[Group_ROOT_DIR,filesep,'Speed_vs_VimMean.fig']);
% saveas(h30,[Group_ROOT_DIR,filesep,'Speed_vs_VimMean.tif']);


h11 = figure(11);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_tracked_branch_d_frames)
    %     text(Group_Pool_tracked_branch_d_frames(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
    plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Tracked Branchs per Frame','fontsize',13);
ylabel('Cell Vim Total Int','fontsize',13);
title({['Tracked Branches per Frame vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')'],...
    ['Correlation: ',num2str(corr(Group_Pool_tracked_branch_d_frames',Group_Pool_whole_cell_vim_totalamount_mean'))]},'fontsize',13);
set(gca,'fontsize',13);

saveas(h11,[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.fig']);
saveas(h11,[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.tif']);
print(h11,'-depsc',[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.eps']);


h12 = figure(12);
plot(Group_Pool_whole_cell_vif_mean_intensity_pat, Group_Pool_branch_duration_array, '.');
ylabel('Branch Duration','fontsize',13);
xlabel('Cell Vim mean Int','fontsize',13);
title('Cell Vim vs Branch Duration','fontsize',13);
set(gca,'fontsize',13);

saveas(h12,[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.fig']);
saveas(h12,[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.tif']);
print(h12,'-depsc',[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.eps']);


h13 = figure(13);
plot(Group_Pool_branch_duration_array(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)), ...
    Group_Pool_branch_vif_mean_intensity(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)),'.');
xlabel('Branch Duration','fontsize',13);
ylabel('Branch Vim mean Int','fontsize',13);
title({['Branch: Duration vs Vim, for branches larger than ',num2str(T_branchsize),' pixels,'],['duration longer than ',num2str(T_branchduration),' frames']},'fontsize',13);
set(gca,'fontsize',13);
saveas(h13,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.fig']);
saveas(h13,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.tif']);
print(h13,'-depsc',[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.eps']);

save([Group_ROOT_DIR,filesep,'branch_analysis_group_results.mat']);


%%
branch_analysis_plotting_VimNms;
branch_analysis_plotting_VimTotalNms;
branch_analysis_plotting_VimFilaDen;
branch_analysis_plotting_VimFilaTotal;


