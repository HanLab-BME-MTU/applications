% function Z = network_analysis_well_orgnized_KS_stat(ML,control_index,feature_index)
% Given the ML(movieList object), load the movies and calculate the KS statistics

% Input: ML: the movieList object for the plate
%        control_index: sepecify which wells are controls
%        feature_index: which features to include

% output: the Z index for each well


load('/project/bioinformatics/Danuser_lab/vimscreen/analysis/tony_zhang/Results/P053_process/row_col_movieList/movieList.mat');

Group_ROOT_DIR='/project/bioinformatics/Danuser_lab/vimscreen/analysis/lding/fromTony/P051052053_process_201601/P053';



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
                try
                    Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
                    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
                    catch
                        col_num = mod(iMD,24);
                        
                        row_num = floor(iMD/24)+1;
                        
                        if(col_num ==0)
                            row_num = row_num-1;
                            col_num =24;
                        end
                        row_col_ID = [char(row_num+64) '_' num2str(col_num,'%02d')];
                        
                        if iChannel ==1
                            
                            for iF = 1 :3
                                filename_short_strs{iF} = ['vim_rrccs_',num2str(row_num,'%02d'),num2str(col_num,'%02d'),num2str(iF,'%01d')];
                            end
                            
                        else
                            for iF = 1 :3
                                filename_short_strs{iF} = ['mt_rrccs_',num2str(row_num,'%02d'),num2str(col_num,'%02d'),num2str(iF,'%01d')];
                            end
                        end
                end
                    
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
                'NA_feature_whole_thisMD');
            
            vim_mean_st_all_cell{1,iMD} = NA_feature_whole_thisMD{1}.mean_st_per_filament_pool;
            vim_sum_st_all_cell{1,iMD} = NA_feature_whole_thisMD{1}.st_per_filament_pool;            
            vim_length_all_cell{1,iMD} = NA_feature_whole_thisMD{1}.length_per_filament_pool;            
            mt_length_all_cell{1,iMD} = NA_feature_whole_thisMD{2}.length_per_filament_pool;
                 
            
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


vim_mean_st_KS=nan(384,384);
vim_sum_st_KS=nan(384,384);
vim_length_KS=nan(384,384);
mt_length_KS=nan(384,384);


for i = 1 :384
    for j = 1:384
        [i j]
        [h,p,ksvalue]=kstest2(vim_mean_st_all_cell{1,i},vim_mean_st_all_cell{1,j},0.05,1);
        vim_mean_st_KS(i,j) = ksvalue;
        [h,p,ksvalue]=kstest2(vim_sum_st_all_cell{1,i},vim_sum_st_all_cell{1,j},0.05,1);
        vim_sum_st_KS(i,j) = ksvalue;
        [h,p,ksvalue]=kstest2(vim_length_all_cell{1,i},vim_length_all_cell{1,j},0.05,1);
        vim_length_KS(i,j) = ksvalue;
        [h,p,ksvalue]=kstest2(mt_length_all_cell{1,i},mt_length_all_cell{1,j},0.05,1);
        mt_length_KS(i,j) = ksvalue;
    end
end

save([Group_ROOT_DIR,filesep,'movieList_plate_4_allKS_gathered.mat'],...
    'vim_mean_st_KS','vim_sum_st_KS','vim_length_KS','mt_length_KS');


vim_mean_st_KS_pn = KS_NP(vim_mean_st_KS);
vim_sum_st_KS_pn = KS_NP(vim_sum_st_KS);
vim_length_KS_pn = KS_NP(vim_length_KS);
mt_length_KS_pn = KS_NP(mt_length_KS);


control_index = C_IDs;
control_pair_numbers = 300;


control_control_index = [];
for iRow = control_index
    for iCol = control_index
        if(iRow~=iCol)
            control_control_index = [control_control_index sub2ind([384,384],iRow,iCol)];
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

control_N_index = [];
for iRow = control_index
    for iCol = N_IDs
        if(iRow~=iCol)
            control_N_index = [control_N_index sub2ind([384,384],iRow,iCol)];
        end
    end
end



A_control_index = [];
for iRow = A_IDs
    for iCol = control_index
        if(iRow~=iCol)
            A_control_index = [A_control_index sub2ind([384,384],iRow,iCol)];
        end
    end
end

control_A_index = [];
for iRow = A_IDs
    for iCol = control_index
        if(iRow~=iCol)
            control_A_index = [control_A_index sub2ind([384,384],iRow,iCol)];
        end
    end
end

m = bootstrp(200, @std, vim_mean_st_KS(control_control_index));    
vim_mean_st_KS_ctrl_bootstrap_std = nanmean(m);
m = bootstrp(200, @std, vim_sum_st_KS(control_control_index));    
vim_sum_st_KS_ctrl_bootstrap_std = nanmean(m);
m = bootstrp(200, @std, vim_length_KS(control_control_index));    
vim_length_KS_ctrl_bootstrap_std = nanmean(m);
m = bootstrp(200, @std, mt_length_KS(control_control_index));    
mt_length_KS_ctrl_bootstrap_std = nanmean(m);

vim_mean_st_Z = vim_mean_st_KS./vim_mean_st_KS_ctrl_bootstrap_std;
vim_sum_st_Z = vim_sum_st_KS./vim_sum_st_KS_ctrl_bootstrap_std;
vim_length_Z = vim_length_KS./vim_length_KS_ctrl_bootstrap_std;
mt_length_Z = mt_length_KS./mt_length_KS_ctrl_bootstrap_std;

m = bootstrp(200, @std, vim_mean_st_KS_pn(control_control_index));    
vim_mean_st_pn_KS_ctrl_bootstrap_std = nanmean(m);
m = bootstrp(200, @std, vim_sum_st_KS_pn(control_control_index));    
vim_sum_st_pn_KS_ctrl_bootstrap_std = nanmean(m);
m = bootstrp(200, @std, vim_length_KS_pn(control_control_index));    
vim_length_pn_KS_ctrl_bootstrap_std = nanmean(m);
m = bootstrp(200, @std, mt_length_KS_pn(control_control_index));    
mt_length_pn_KS_ctrl_bootstrap_std = nanmean(m);

vim_mean_st_pn_Z = vim_mean_st_KS_pn./vim_mean_st_pn_KS_ctrl_bootstrap_std;
vim_sum_st_pn_Z = vim_sum_st_KS_pn./vim_sum_st_pn_KS_ctrl_bootstrap_std;
vim_length_pn_Z = vim_length_KS_pn./vim_length_pn_KS_ctrl_bootstrap_std;
mt_length_pn_Z = mt_length_KS_pn./mt_length_pn_KS_ctrl_bootstrap_std;

Z = {vim_mean_st_Z,vim_sum_st_Z,vim_length_Z,mt_length_Z,vim_mean_st_pn_Z,vim_sum_st_pn_Z,vim_length_pn_Z,mt_length_pn_Z};


for iRow = 1:384

m = bootstrp(200, @mean, vim_mean_st_Z(iRow,control_index));    
vim_mean_st_Z_bootstrap_mean(iRow) = nanmean(m);
m = bootstrp(200, @std, vim_sum_st_Z(iRow,control_index));   
vim_sum_st_Z_bootstrap_mean(iRow) = nanmean(m);
m = bootstrp(200, @std, vim_length_Z(iRow,control_index));     
vim_length_Z_bootstrap_mean(iRow) = nanmean(m);
m = bootstrp(200, @std, mt_length_Z(iRow,control_index));   
mt_length_Z_bootstrap_mean(iRow) = nanmean(m);

end
  

for iRow = 1:384

m = bootstrp(20, @mean, vim_mean_st_pn_Z(iRow,control_index));    
vim_mean_st_pn_Z_bootstrap_mean(iRow) = nanmean(m);
m = bootstrp(20, @std, vim_sum_st_pn_Z(iRow,control_index));   
vim_sum_st_pn_Z_bootstrap_mean(iRow) = nanmean(m);
m = bootstrp(20, @std, vim_length_pn_Z(iRow,control_index));     
vim_length_pn_Z_bootstrap_mean(iRow) = nanmean(m);
m = bootstrp(20, @std, mt_length_pn_Z(iRow,control_index));   
mt_length_pn_Z_bootstrap_mean(iRow) = nanmean(m);

end


CCCC1=nan(24,16);
CCCC1(:)=vim_mean_st_pn_Z_bootstrap_mean(:);
% CCCC(C_IDs)=0;
figure;
plot_color_dot_from_matrix(flipud(CCCC1'),20,-2, 2,1)

CCCC2=nan(24,16);
CCCC2(:)=vim_sum_st_pn_Z_bootstrap_mean(:);
% CCCC(C_IDs)=0;
figure;
plot_color_dot_from_matrix(flipud(CCCC2'),20,-2, 2,1)

CCCC3=nan(24,16);
CCCC3(:)=vim_length_pn_Z_bootstrap_mean(:);
% CCCC(C_IDs)=0;
figure;plot_color_dot_from_matrix(flipud(CCCC3'),20,-2, 2,1)

CCCC4=nan(24,16);
CCCC4(:)=mt_length_pn_Z_bootstrap_mean(:);
% CCCC(C_IDs)=0;
figure;plot_color_dot_from_matrix(flipud(CCCC4'),20,-2, 2,1);  set(gca,'XTick',0.5:24-0.5);

figure;plot_color_dot_from_matrix(flipud(abs(CCCC1')+abs(CCCC2')+abs(CCCC3')+abs(CCCC4')),20, 0, 5, 1);


