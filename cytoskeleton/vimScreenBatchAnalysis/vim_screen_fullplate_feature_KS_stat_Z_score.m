function [Z_maps, KS_maps] = vim_screen_fullplate_feature_KS_stat_Z_score(plate_analysis_folder, feature_index, control_index, positive_index,feature_gathered_folder, ks_results_saving_folder)
% Function for calculate the KS statistics for the network analysis resulting features for the full plate of 384 wells.

% Liya Ding
% Matlab 2015a
% 2016.02

% note: this version will take the default setting and default index, it
% might or might not work with different setting and index, need to
% reconfirm later

% Note: due to the size of the data and analysis, the raw images and some
% results will be removed after the FilamentAnalysisPackage and network
% analysis job finishes. So the movieData structure is broken. To deal with
% this, given the vim screen struture is fixed, the folders here are hard coded.
% For the future, if the structure updates, the code need to be updated as well.

% Input: plate_analysis_folder :              The root folder of the plate analysis. (String of the folder name).
%        feature_index:                       The features interested. If not given, use all 28 possible features.
%                                                ( note that after re-organization, there are 33 features.)
%        control_index:                       Sepecify which wells are negative controls;
%                                                currently, if not given, use the control plates index
%        positive_index:                      Sepecify which wells are positive controls;
%                                                currently, if not given, use the control plates index
%        feature_gathered_folder:             The folder where the plate grouped feature mat files are saved
%                                                If not given, this will be the movieList folder.
%        ks_results_saving_folder:            The folder to save all the results for KS tests results.
%                                                If not given, this will be the movieList folder.

% output: Z_maps:                                  includes a few fields
%                                           The z score maps per channel for feature 
%                                                    Z_maps.feature_Z_bootstrap_mean_cell
%                                                    Z_maps.feature_bootstrap_std_cell
%                                           the z score when combined and whether they are hits based on the thresholds
%                                                    Z_select_abs_thisplate
%                                                    Z_maps.Z_L1
%                                                    Z_maps.Z_L2
%                                                    Z_maps.Z_hit_L1
%                                                    Z_maps.Z_hit_L1_lowerbound
%                                                    Z_maps.Z_hit_L2
%                                                    Z_maps.Z_hit_L2_lowerbound
%   KS_maps:
%               KS_maps.feature_PN_KS_cell:   the ks value per well pair, positive and negative with dominant value
%               KS_maps.feature_KS_cell:   the ks value per well pair, original

% initialization
KS_maps = [];
Z_maps = [];

% if no feature index is given, use all features
if(nargin<2)
    feature_index=ones(28,1);
end

% if no negative control well index given, use default
if(nargin<3)
    A_IDs = [24*(1:14)+1 3*24+10];
    N_IDs = [24*(1:14)+24 9*24+15];
    C_IDs = 1:384;
    C_IDs = setdiff(setdiff(C_IDs,A_IDs),N_IDs);
    control_index = C_IDs;    
end

% if no positive control well index given, use default
if(nargin<4)
    A_IDs = [24*(1:14)+1 3*24+10];
    N_IDs = [24*(1:14)+24 9*24+15];
    C_IDs = 1:384;
    C_IDs = setdiff(setdiff(C_IDs,A_IDs),N_IDs);
    positive_index = A_IDs;
end

% if no output saving folder is given, use the movieList folder under plate analysis folder
if(nargin<5)
    feature_gathered_folder = [plate_analysis_folder, filesep,'row_col_movieList'];
end

% if no output saving folder is given, use the movieList folder under plate analysis folder
if(nargin<6)
    ks_results_saving_folder = [plate_analysis_folder, filesep,'row_col_movieList'];
end

movieNumber = 384; % fixed

ML_name = [plate_analysis_folder, filesep,'row_col_movieList',filesep,'movieList.mat'];

% load directly, not movieList.load since the list of MD are broken.
load(ML_name); % then we have ML

movieNumber = numel(ML.movieDataFile_);

for iMD  = 1 : movieNumber
    
    % otherwise load one by one
    NA_feature_thisMD = cell(2,1);
    Identifier_thisMD = cell(2,1);
    NA_feature_whole_thisMD = cell(2,1);
    
    try
        % load this movie
        load(ML.movieDataFile_{iMD});
    catch
        disp('The MD file is missing');
        continue;
    end
    
    % Now MD in workspace
    %check index
    display_msg_flag=0;
    package_process_ind_script;
    
    nChannel = numel(MD.channels_);
    nFrame = MD.nFrames_;
    
    for iChannel = 1 : 2
        display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);
        
        outdir = [MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(iChannel),filesep,'analysis_results'];
        Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
        filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
        
        % check if this folder exist
        if(exist(outdir,'dir'))
            % if it exist, try to do the branch analysis
            for iFrame = 1 : nFrame
                display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
                
                %
                Identifier_thisMD{iChannel, iFrame} = filename_short_strs{iFrame};
                
            end% end of a frame
            
        end % end of if the analysis for output folder exist
        
    end % end of a Channel
    
    load([feature_gathered_folder,filesep,Identifier_thisMD{iChannel, iFrame},'_movieData_well_NA_results_pool_gathered.mat'],...
        'Identifier_thisMD',...
        'ChMP_feature_thisMD',...
        'NA_feature_whole_thisMD',...
        'valid_feature_index',...
        'valid_requested_feature_index');
    
    for iChannel = 1 : 2
        for iF = 1:33
            
            if(valid_requested_feature_index(iF)>0)
                switch iF
                    case 1
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.straightness_per_filament_pool;
                    case 2
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.length_per_filament_pool;
                    case 3
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.pixel_number_per_filament_pool;
                    case 4
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.density_filament;
                    case 5
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.scrabled_density_filament;
                    case 6
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display;
                    case 7
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.orientation_pixel_pool_display_center;
                    case 8
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.intensity_per_filament_pool;
                    case 9
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.mean_intensity_per_filament_pool;
                    case 10
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.intensity_per_fat_filament_pool;
                    case 11
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.mean_intensity_per_fat_filament_pool;
                    case 12
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.scale_per_filament_pool;
                    case 13
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.st_per_filament_pool;
                    case 14
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.mean_st_per_filament_pool;
                    case 15
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.st_per_fat_filament_pool;
                    case 16
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.mean_st_per_fat_filament_pool;
                    case 17
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.filament_mean_curvature;
                        end
                    case 18
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.curvature_per_pixel_pool;
                        end
                    case 19
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool;
                        end
                    case 20
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool;
                        end
                    case 21
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool;
                        end
                    case 22
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool;
                        end
                    case 23
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool;
                        end
                    case 24
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool;
                        end
                    case 25
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.nmssum_ratio_pool_autoDpc;
                        end
                    case 26
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.nmsmean_ratio_pool_autoDpc;
                        end
                    case 27
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.intsum_ratio_pool_autoDpc;
                        end
                    case 28
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.intmean_ratio_pool_autoDpc;
                        end
                    case 29
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.filasum_ratio_pool_autoDpc;
                        end
                    case 30
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.filamean_ratio_pool_autoDpc;
                        end
                    case 31
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.Centripetal_fila;
                        end
                    case 32
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.Centripetal_pixel;
                        end
                    case 33
                        try
                            feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.number_of_nucleus;
                        end
                    otherwise
                        feature_pool_all_cell{iChannel,iMD}{iF}= NA_feature_whole_thisMD{iChannel}.straightness_per_filament_pool;
                end
            end
        end
    end
end% end of a MD

%%
A_IDs = [24*(1:14)+1 3*24+10];
N_IDs = [24*(1:14)+24 9*24+15];
C_IDs = 1:384;
C_IDs = setdiff(setdiff(C_IDs,A_IDs),N_IDs);

feature_KS_cell = cell(2,33);
feature_PN_KS_cell = cell(2,33);

for iF = 1 : 33
    if(valid_requested_feature_index(iF)>0)
        feature_KS_cell{iChannel,iF} = nan(384,384);
        for iChannel = 1 :2
            for i = 1 :384
                for j = 1:384
                    [i j]
                    [h,p,ksvalue]=kstest2(feature_pool_all_cell{iChannel,i}{iF},feature_pool_all_cell{iChannel,j}{iF},0.05,1);
                    feature_KS_cell{iChannel,iF}(i,j) = ksvalue;
                end
            end
            feature_PN_KS_cell{iChannel,iF} = KS_NP(feature_KS_cell{iChannel,iF});
        end
    end
end

for iF = 1 : 33
    if(valid_requested_feature_index(iF)>0)
        for iChannel = 1 :2
            feature_PN_KS_cell{iChannel,iF} = KS_NP(feature_KS_cell{iChannel,iF});
        end
    end
end


save([ks_results_saving_folder,filesep,'movieList_plate_allKS_gathered.mat'],...
    'feature_KS_cell', 'feature_PN_KS_cell');

KS_maps.feature_PN_KS_cell = feature_PN_KS_cell;
KS_maps.feature_KS_cell = feature_KS_cell;


%% the pairs index
control_pair_numbers = 300;

% get the control pair indices
control_control_index = [];
for iRow = control_index
    for iCol = control_index
        if(iRow~=iCol)
            control_control_index = [control_control_index sub2ind([384,384],iRow,iCol)];
        end
    end
end


% get the control pair indices
A_control_index = [];
for iRow = A_IDs
    for iCol = control_index
        if(iRow~=iCol)
            A_control_index = [A_control_index sub2ind([384,384],iRow,iCol)];
        end
    end
end
for iRow = control_index
    for iCol = A_IDs
        if(iRow~=iCol)
            A_control_index = [A_control_index sub2ind([384,384],iRow,iCol)];
        end
    end
end

N_control_index = [];
for iRow = N_IDs
    for iCol = control_index
        if(iRow~=iCol)
            N_control_index = [N_control_index sub2ind([384,384],iRow,iCol)];
        end
    end
end
for iRow = control_index
    for iCol = N_IDs
        if(iRow~=iCol)
            N_control_index = [N_control_index sub2ind([384,384],iRow,iCol)];
        end
    end
end

%% calculate the Z scores

% Z = KS/bootstrap_std

sum_Z = zeros(16,24);
feature_Z_cell=cell(2,33);
feature_Z_bootstrap_mean_cell=cell(2,33);
feature_bootstrap_std_cell=cell(2,33);

for iChannel = 1 :2
    for iF = 1 : 33
        if(valid_requested_feature_index(iF)>0)
            
            % get the std with bootstrap from the control pairs
            m = bootstrp(200, @std, feature_PN_KS_cell{iChannel,iF}(control_control_index));
            feature_bootstrap_std = nanmean(m);
            feature_bootstrap_std_cell{iChannel,iF} = feature_bootstrap_std;
            
            feature_Z_cell{iChannel,iF} = feature_PN_KS_cell{iChannel,iF}./feature_bootstrap_std;
            
            for iRow = 1:384
                m = bootstrp(20, @mean, feature_Z_cell{iChannel,iF}(iRow,control_index));
                feature_Z_bootstrap_mean_cell{iChannel,iF}(iRow) = nanmean(m);
            end
            
            CCCC1=nan(24,16);
            CCCC1(:)=feature_Z_bootstrap_mean_cell{iChannel,iF}(:);
            figure;
            plot_color_dot_from_matrix(flipud(CCCC1'),20,-2, 2,1);
            close;
            
            sum_Z = sum_Z + abs(CCCC1');
        end
    end
end

figure;plot_color_dot_from_matrix(flipud(sum_Z),20, 0, 45, 1);
close;

%% Z score selection, combination and gating
% these two are the slected features for now
selected_index_F =[
     1
     8
     9
    10
    11
    13
    15
     1
     8
    10
    13
    15];
    
selected_index_C = [ 
     1
     1
     1
     1
     1
     1
     1
     2
     2
     2
     2
     2];
 
 
Z_select_thisplate = [];
for ind = 1 :12
    iC = selected_index_C(ind);
    iF = selected_index_F(ind);
    Z_select_thisplate =[Z_select_thisplate; feature_Z_bootstrap_mean_cell{iC,iF};];    
end

save([ks_results_saving_folder,filesep,'movieList_plate_Z_results.mat'],...
    'Z_select_thisplate', 'feature_Z_bootstrap_mean_cell','selected_index_C','selected_index_F');

Z_select_norm_thisplate = sqrt(sum(Z_select_thisplate.^2,1));
Z_select_sum_thisplate  = (sum(abs(Z_select_thisplate),1));
Z_select_abs_thisplate = abs(Z_select_thisplate);
Z_select_min_thisplate = min(abs(Z_select_thisplate));


CCCC1=nan(24,16);
CCCC1(:) = Z_select_norm_thisplate(:);
Z_hit_L2 = flipud(CCCC1')>sqrt(2.5*2.5*12);


% add additional lowerbound for each Z value, set here as 1.25
CCCC1=nan(24,16);
CCCC1(:) = Z_select_min_thisplate(:);
Z_hit_lowerbound = flipud(CCCC1')>1.25;

Z_hit_L2_lowerbound = Z_hit_L2.*Z_hit_lowerbound; 

h = figure;
plot_color_dot_from_matrix(Z_hit_L2,20,0,1,1);
title('Plate: Z score L^2 norm gated','fontsize',13);
  set(gca,'fontsize',13);
  saveas(h,[ks_results_saving_folder,filesep,'Z_hit_L2.tif']);
     
 h = figure;
plot_color_dot_from_matrix(Z_hit_L2_lowerbound,20,0,1,1);
title('Plate: Z score L^2 norm gated with lowerbound on each Z','fontsize',13);
  set(gca,'fontsize',13);
 saveas(h,[ks_results_saving_folder,filesep,'Z_hit_L2_lowerbound.tif']);
     
CCCC1=nan(24,16);
CCCC1(:)=Z_select_sum_thisplate(:);
Z_hit_L1 = flipud(CCCC1')>2.5*12;

% add additional lowerbound for each Z value, set here as 1.25
Z_hit_L1_lowerbound = Z_hit_L1.*Z_hit_lowerbound; 

h = figure;
plot_color_dot_from_matrix(Z_hit_L1,20,0,1,1);
title('Plate: Z score L^1 norm gated','fontsize',13);
set(gca,'fontsize',13);
  
saveas(h,[ks_results_saving_folder,filesep,'Z_hit_L1.tif']);
     

h = figure;
plot_color_dot_from_matrix(Z_hit_L1_lowerbound,20,0,1,1);
title('Plate: Z score L^1 norm gated with lowerbound on each Z','fontsize',13);
set(gca,'fontsize',13);
saveas(h,[ks_results_saving_folder,filesep,'Z_hit_L1_lowerbound.tif']);


Z_maps.Z_L2 = Z_select_norm_thisplate;
Z_maps.Z_L1  = Z_select_sum_thisplate;
Z_maps.Z_select_abs_thisplate = Z_select_abs_thisplate;

Z_maps.Z_hit_L1 = Z_hit_L1;
Z_maps.Z_hit_L1_lowerbound = Z_hit_L1_lowerbound;
  Z_maps.Z_hit_L2 = Z_hit_L2;
  Z_maps.Z_hit_L2_lowerbound = Z_hit_L2_lowerbound;
  Z_maps.feature_Z_bootstrap_mean_cell = feature_Z_bootstrap_mean_cell;
  Z_maps.feature_bootstrap_std_cell = feature_bootstrap_std_cell;
  
    save([ks_results_saving_folder,filesep,'movieList_plate_Z_cells.mat']);

