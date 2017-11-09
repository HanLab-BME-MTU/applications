% 
% ML_name_cell=[];
% ML_name_cell{1} = '/project/cellbiology/gdanuser/vimentin/tonyzhang/screen_data/Nikon/data_analysis/P042_process/row_col_movieList/movieList.mat';
% 
% feature_index=zeros(28,1);
% feature_index(1:3)=1;
% feature_index(8:16)=1;
% 
% Group_ROOT_DIR='/project/cellbiology/gdanuser/vimentin/ding/fromTony/P042_process_20151115';
% 
% 
% close all;
% 
% nList = length(ML_name_cell);
% 
% vim_mean_st_all_cell = cell(1,384);
% vim_sum_st_all_cell = cell(1,384);
% vim_length_all_cell = cell(1,384);
% 
% mt_length_all_cell = cell(1,384);
% 
% for iML = 1 : nList
%     
%    try
%        load(ML_name_cell{iML});
%    catch
%        display('Error in ML loading.');
%        break;
%    end
%    
%    ML_ROOT_DIR = Group_ROOT_DIR;
%     
%     % the number of movies
%     movieNumber =  length(ML.movieDataFile_);
%     
%     NA_feature_thisML = cell(2,1);
%     Identifier_thisML = cell(2,1);
%     NA_feature_whole_thisMD = cell(2,1);
%     
%     % if the batch results exist, load it
%     if 0
%         
%         % this part commented because the result for a whole ML will be too
%         % large to store
% %         (exist([ML_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'], 'file'))
% %         load([ML_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'],...
% %             'NA_feature_thisML','Identifier_thisML');
%     else
%         for iMD  = 1 : movieNumber
%             
%                       
%             % otherwise load one by one
%             NA_feature_thisMD = cell(2,1);
%             Identifier_thisMD = cell(2,1);
%             NA_feature_whole_thisMD = cell(2,1);
%   
%             if 0
% %                 (exist([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'], 'file'))
% %                 load([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'],...
% %                     'NA_feature_thisMD','Identifier_thisMD');
%             else                
%                 try
%                     % load this movie
%                     load(ML.movieDataFile_{iMD});
%                 catch
%                     disp('The MD file is missing');
%                     continue;
%                 end
%                 
%                 % Now MD in workspace
%                 %check index
%                 display_msg_flag=0;
%                 package_process_ind_script;
%                 
%                 MD_ROOT_DIR = MD.outputDirectory_;
%                 
%                 nChannel = numel(MD.channels_);
%                 nFrame = MD.nFrames_;
%                 
%                 for iChannel = 1 : 2
%                     display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);
% %                     outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
%                     
%                     outdir = [MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(iChannel),filesep,'analysis_results'];
%                     Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
%                     filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
%                     
%                     % check if this folder exist
%                     if(exist(outdir,'dir'))
%                         % if it exist, try to do the branch analysis
%                         %                         try
%                         
%                         for iFrame = 1 : nFrame
%                             display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
% 
%                             %
%                             Identifier_thisMD{iChannel, iFrame} = filename_short_strs{iFrame};      
%                                                      
%                         end% end of a frame
%                         
%                     end % end of if the analysis for output folder exist
%                     
%                                                        
%                 end % end of a Channel
%                
%                 Identifier_thisMD{iChannel, iFrame}
%                 
%                 load([Group_ROOT_DIR,filesep,Identifier_thisMD{iChannel, iFrame},'_movieData_well_NA_results_pool_gathered.mat'],...
%                     'Identifier_thisMD',...
%                     'ChMP_feature_thisMD',...
%                     'NA_feature_whole_thisMD');
%                 
%                 
%                 
%                 vim_mean_st_all_cell{1,iMD} = NA_feature_whole_thisMD{1}.mean_st_per_filament_pool;
%                 vim_sum_st_all_cell{1,iMD} = NA_feature_whole_thisMD{1}.st_per_filament_pool;
%                                
%                 vim_length_all_cell{1,iMD} = NA_feature_whole_thisMD{1}.length_per_filament_pool;
%                 
%                 mt_length_all_cell{1,iMD} = NA_feature_whole_thisMD{2}.length_per_filament_pool;
%            end % end of if previous gathering exists for MD
%             
%          
%         end  % end of a MD        
%         
%         
%     end     % end of if previous gathering exists FOR a ML 
%   
% end % end of a ML
% 
%   %if mod(iMD,24)==0
%      save([Group_ROOT_DIR,filesep,'movieList_plate_4_NA_results_gathered.mat'],...
%                     'vim_mean_st_all_cell', 'vim_sum_st_all_cell','vim_length_all_cell','mt_length_all_cell');
%     
%    %        end

 A_IDs = [24*(1:14)+1 3*24+10];
 N_IDs = [24*(1:14)+24 9*24+15];
 C_IDs = 1:384;
 C_IDs = setdiff(setdiff(C_IDs,A_IDs),N_IDs);
%  
%  
% %  vim_mean_st_KS=nan(384,384);
% %   vim_sum_st_KS=nan(384,384);
% %  vim_length_KS=nan(384,384);
% %   mt_length_KS=nan(384,384);
% %  
%  
% 
%   for i = 1 :383
%       
%     for j = i+1:384
%         [i j]
%         [h,p,ksvalue]=kstest2(vim_mean_st_all_cell{1,i},vim_mean_st_all_cell{1,j});
%         vim_mean_st_KS(i,j) = ksvalue;
%           [h,p,ksvalue]=kstest2(vim_sum_st_all_cell{1,i},vim_sum_st_all_cell{1,j});
%         vim_sum_st_KS(i,j) = ksvalue;
%            [h,p,ksvalue]=kstest2(vim_length_all_cell{1,i},vim_length_all_cell{1,j});
%         vim_length_KS(i,j) = ksvalue;
%         [h,p,ksvalue]=kstest2(mt_length_all_cell{1,i},mt_length_all_cell{1,j});
%         mt_length_KS(i,j) = ksvalue;
%       
%     end
%      if mod(i,12)==0
%     save([Group_ROOT_DIR,filesep,'movieList_plate_4_KS_gathered_1121.mat'],...
%                    'vim_mean_st_KS','vim_sum_st_KS','vim_length_KS','mt_length_KS');
%      end
% end
% 
%   save([Group_ROOT_DIR,filesep,'movieList_plate_4_KS_gathered_1121.mat'],...
%                    'vim_mean_st_KS','vim_sum_st_KS','vim_length_KS','mt_length_KS');
        Group_ROOT_DIR = '/project/cellbiology/gdanuser/vimentin/ding/fromTony';
     mt_length_control_big_pool =[];
      vim_length_control_big_pool =[];
     vim_mean_st_control_big_pool =[];    
      vim_sum_st_control_big_pool =[];
              
   for i = 1 :384
       
        mt_length_control_big_pool =[mt_length_control_big_pool; mt_length_all_cell{1,i}];
      vim_length_control_big_pool =[vim_length_control_big_pool; vim_length_all_cell{1,i}];
     vim_mean_st_control_big_pool =[vim_mean_st_control_big_pool; vim_mean_st_all_cell{1,i}];    
      vim_sum_st_control_big_pool =[vim_sum_st_control_big_pool; vim_sum_st_all_cell{1,i}];
     
     if mod(i,24)==0 
         save([Group_ROOT_DIR,filesep,'movieList_plate_4_controlbigpool.mat'],...
              'mt_length_control_big_pool',      'vim_length_control_big_pool',     'vim_mean_st_control_big_pool',          'vim_sum_st_control_big_pool');
  
     end
   end             
               
    save([Group_ROOT_DIR,filesep,'movieList_plate_4_controlbigpool.mat'],...
              'mt_length_control_big_pool',      'vim_length_control_big_pool',     'vim_mean_st_control_big_pool',          'vim_sum_st_control_big_pool');
              
   
    for j = 1 :384
        
        j
        [h,p,ksvalue]=kstest2(vim_mean_st_control_big_pool,vim_mean_st_all_cell{1,j});
        vim_mean_st_KS_toCtrl(1,j) = ksvalue;
        [h,p,ksvalue]=kstest2(vim_sum_st_control_big_pool,vim_sum_st_all_cell{1,j});
        vim_sum_st_KS_toCtrl(1,j) = ksvalue;
        [h,p,ksvalue]=kstest2(vim_length_control_big_pool,vim_length_all_cell{1,j});
        vim_length_KS_toCtrl(1,j) = ksvalue;
        [h,p,ksvalue]=kstest2(mt_length_control_big_pool,mt_length_all_cell{1,j});
        mt_length_KS_toCtrl(1,j) = ksvalue;
       
      
     if mod(j,24)==0 
         save([Group_ROOT_DIR,filesep,'movieList_plate_4_KS_tocontrolbigpool.mat'],...
              'vim_mean_st_KS_toCtrl',      'vim_sum_st_KS_toCtrl', ...
              'vim_length_KS_toCtrl',          'mt_length_KS_toCtrl');
  
     end
   end       
          
          