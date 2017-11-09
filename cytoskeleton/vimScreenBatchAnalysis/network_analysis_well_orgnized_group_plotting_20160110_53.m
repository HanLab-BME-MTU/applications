
ML_name_cell=[];
ML_name_cell{1} = '/project/bioinformatics/Danuser_lab/vimscreen/analysis/tony_zhang/Results/P053_process/row_col_movieList/movieList.mat';

feature_index=zeros(28,1);
feature_index(1:28)=1;

Group_ROOT_DIR='/project/bioinformatics/Danuser_lab/vimscreen/analysis/lding/fromTony/P051052053_process_201601/P053';


close all;

nList = length(ML_name_cell);

for iML = 1 : nList
    
   try
       load(ML_name_cell{iML});
   catch
       display('Error in ML loading.');
       break;
   end
   
   ML_ROOT_DIR = Group_ROOT_DIR;
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    NA_feature_thisML = cell(2,1);
    Identifier_thisML = cell(2,1);
    NA_feature_whole_thisMD = cell(2,1);
    
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
                
                MD_ROOT_DIR = MD.outputDirectory_;
                
                nChannel = numel(MD.channels_);
                nFrame = MD.nFrames_;
                
                for iChannel = 1 : 2
                    display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);
%                     outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
                    
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
                        %                         try
                        
                        for iFrame = 1 : nFrame
                            display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
                            load_ch_frame_flag =0;
                            
                            if(exist([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),...
                                    '_frame_',num2str(iFrame),'.mat'],'file'))
                                load([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),...
                                    '_frame_',num2str(iFrame),'.mat'],'output_feature');
                                load_ch_frame_flag=1;
                            end
                            
                            if(exist([outdir,filesep,'network_feature_ch_',num2str(iChannel),'_',...
                                    filename_short_strs{iFrame},'.mat'],'file'))
                                load([outdir,filesep,'network_feature_ch_',num2str(iChannel),'_',...
                                    filename_short_strs{iFrame},'.mat'], 'output_feature');
                                load_ch_frame_flag = 1;
                            end
                            
                            if(load_ch_frame_flag==0)
                                continue;
                            end
                             
                            %arrange the features, all into column vectors(density, keep only ROI part)                            
                            output_feature_arranged = arrange_network_features(output_feature);
                             
                            % organize profiles
                            output_feature_reorganized = organize_profiles(output_feature_arranged,feature_index);

                            % update the index as well
                            feature_index_reorganized = feature_index;
                            
                            if(feature_index(19) == 1 && ~isempty(output_feature_reorganized.profileAllCell))
    
                                feature_index_reorganized(19:30) = 1; 
                                feature_index_reorganized(31:33) = feature_index(21:23);                                 
                            else
                                feature_index_reorganized(19:30) = 0; 
                                feature_index_reorganized(31:33) = feature_index(21:23); 
                            end                   
                            
                             % find the feature index that is valid( with valid calculated values)
                          valid_feature_index = find_valid_feature_index(output_feature_reorganized,ones(33,1));
                             % find the feature index that is valid( with
                            % valid calculated values) and requested
                            valid_requested_feature_index = valid_feature_index.*feature_index_reorganized;
                           
                            % keep this into cells
                            NA_feature_thisMD{iChannel, iFrame} = output_feature_reorganized;
                            
                            % 
                            Identifier_thisMD{iChannel, iFrame} = filename_short_strs{iFrame};      
                            
                            % get the movie pools 
                             % not using function, the variables are way too large in size
%                            NA_feature_whole_thisMD = put_pool_together(NA_feature_whole_thisMD, NA_feature_thisMD,iChannel);                 

                          put_frame_pool_into_MD_pool;
                            
                                                      
                        end% end of a frame
                        
                    end % end of if the analysis for output folder exist
                    
                    
                    % since each movie is each well, put the pool together, and get the well percentiles, (for a channel)
                    tic
                    ChMP_feature_thisMD{iChannel} = extract_mean_percentiles_well(NA_feature_whole_thisMD, iChannel, valid_requested_feature_index);
                    toc
                                        
                                                                           
                end % end of a Channel
               
                
                save([Group_ROOT_DIR,filesep,Identifier_thisMD{iChannel, iFrame},'_movieData_well_NA_results_pool_gathered.mat'],...
                    'Identifier_thisMD',...
                    'ChMP_feature_thisMD',...
                    'NA_feature_whole_thisMD',...
                    'valid_feature_index',...
                    'valid_requested_feature_index');

            end % end of if previous gathering exists for MD
            
            Identifier_thisML{1,iMD} = Identifier_thisMD;
            ChMP_feature_thisML{1,iMD} = ChMP_feature_thisMD;

        end  % end of a MD        
        
         save([Group_ROOT_DIR,filesep,'movieList_plate_NA_results_gathered.mat'],...
                    'Identifier_thisML',...
                    'ChMP_feature_thisML');
    end     % end of if previous gathering exists FOR a ML    
    
    
%     %% network analysis plotting for vim screen
%     plate_network_feature_plotting(Identifier_thisML,...
%                     CFMP_feature_ordered_thisML,...
%                     CFMP_feature_thisML,...
%                     ChMP_feature_thisML,ChMP_feature_perwell_thisML);
%     
    
    
end % end of a ML

