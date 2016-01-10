% function Z = network_analysis_well_orgnized_KS_stat(ML,control_index,feature_index)
% Given the ML(movieList object), load the movies and calculate the KS statistics

% Input: ML: the movieList object for the plate
%        control_index: sepecify which wells are controls
%        feature_index: which features to include

% output: the Z index for each well


load('/project/bioinformatics/Danuser_lab/vimscreen/analysis/tony_zhang/Results/P052_process/row_col_movieList/movieList.mat');

Group_ROOT_DIR='/project/bioinformatics/Danuser_lab/vimscreen/analysis/lding/fromTony/P051052053_process_201601/P052';



vim_mean_st_all_cell = cell(1,384);
vim_sum_st_all_cell = cell(1,384);
vim_length_all_cell = cell(1,384);
mt_length_all_cell = cell(1,384);

% the number of movies
movieNumber =  length(ML.movieDataFile_);


% if the batch results exist, load it
if 0
    
    % this part commented because the result for a whole ML will be too
    % large to store
    %         (exist([ML_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'], 'file'))
    %         load([ML_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'],...
    %             'NA_feature_thisML','Identifier_thisML');
else
    for iMD  = 1 : movieNumber
        
        % otherwise load one by one
        NA_feature_thisMD = cell(2,1);
        Identifier_thisMD = cell(2,1);
        NA_feature_whole_thisMD = cell(2,1);
        
        if 0
            %                 (exist([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'], 'file'))
            %                 load([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'],...
            %                     'NA_feature_thisMD','Identifier_thisMD');
        else
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
            
            load([Group_ROOT_DIR,filesep,Identifier_thisMD{iChannel, iFrame},'_movieData_well_NA_results_pool_gathered.mat'],...
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
        end % end of if previous gathering exists for MD
        
        
    end  % end of a MD
    
    
end     % end of if previous gathering exists FOR a ML


%if mod(iMD,24)==0
save([Group_ROOT_DIR,filesep,'movieList_plate_4_NA_results_gathered.mat'],...
    'vim_mean_st_all_cell', 'vim_sum_st_all_cell','vim_length_all_cell','mt_length_all_cell');

%        end

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
        end
        
        feature_PN_KS_cell{iChannel,iF} = KS_NP(feature_KS_cell{iChannel,iF});
    end
end

save([Group_ROOT_DIR,filesep,'movieList_plate_allKS_gathered.mat'],...
    'feature_KS_cell', 'feature_PN_KS_cell');


%% calculate the Z scores

control_index = C_IDs;
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

% Z = KS/bootstrap_std

sum_Z = zeros(24,16);
feature_Z_cell=cell(2,33);
feature_Z_bootstrap_mean_cell=cell(2,33);

for iChannel = 1 :2
    for iF = 1 : 33
        if(valid_requested_feature_index(iF)>0)
            
            % get the std with bootstrap from the control pairs
            m = bootstrp(200, @std, vim_mean_st_KS_pn(control_control_index));
            feature_bootstrap_std = nanmean(m);
            
            feature_Z_cell{iChannel,iF} = feature_PN_KS_cell{iChannel,iF}./feature_bootstrap_std;
            
            for iRow = 1:384
                m = bootstrp(20, @mean, feature_Z_cell{iChannel,iF}(iRow,control_index));
                feature_Z_bootstrap_mean_cell{iChannel,iF} = nanmean(m);
            end
            
            CCCC1=nan(24,16);
            CCCC1(:)=feature_Z_bootstrap_mean_cell{iChannel,iF}(:);
            figure;
            plot_color_dot_from_matrix(flipud(CCCC1'),20,-2, 2,1);
            
            sum_Z = sum_Z + abs(CCCC1');
        end
    end
end

figure;plot_color_dot_from_matrix(flipud(sum_Z),20, 0, 5, 1);

save([Group_ROOT_DIR,filesep,'movieList_plate_Z_all.mat']);

