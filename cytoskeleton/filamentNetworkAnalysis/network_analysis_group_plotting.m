% function NA_Group = network_analysis_group_plotting(ML_name_cell, feature_index, Group_ROOT_DIR)
% % function to do network feature plotting for a whole movielist and then
% % for a list of ML
% % Liya Ding, March, 2014
% %
% % Input:
% %   ML_array:     The array of movieList objects in this group to be analyzed
% 
% 
% if(nargin<2)
%     feature_index=ones(18,1);
% end
% 
% if(nargin<3)
%     Group_ROOT_DIR=[];
% end
ML_name_cell=[];
ML_name_cell{1} = '/project/cellbiology/gdanuser/vimentin/ding/fromTony/screen_20141204_3column_test/threecolumn_movieList/movieList.mat';
feature_index=ones(18,1);
Group_ROOT_DIR=[];


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
            
            MD=MovieData.load(ML.movieDataFile_{iMD});
             
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
                            
                            if(exist([outdir,filesep,'network_analysis_feature_ch_',num2str(iChannel),...
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
                            
                            NA_feature_thisMD{iChannel, iFrame} = output_feature;
                            Identifier_thisMD{iChannel, iFrame} = filename_short_strs{iFrame};                            
                            CFMP_feature_thisMD{iChannel, iFrame} = extract_mean_percentiles_one(output_feature, feature_index);
                            
                                                        
                        end% end of a frame
                        
                    end % end of if the analysis for output folder exist
                    
                    ChMP_feature_thisMD{iChannel} = extract_mean_percentiles_cell(NA_feature_thisMD,iChannel, feature_index);
                    
                    % remove it after use, since it is too big                    
                    for iFrame = 1 : nFrame
                        NA_feature_thisMD{iChannel, iFrame}=[];
                    end
                    
                    CFMP_feature_ordered_thisMD{iChannel} = nan(18,8,nFrame);
                    ALL_cat = cat(1,CFMP_feature_thisMD{:});
                    
                    for iF = 1 :18
                        for iP = 1 :8
                           CFMP_feature_ordered_thisMD{iChannel}(iF,iP,1:nFrame) = ALL_cat(iF:18:end,iP);
                        end
                    end                 
                end % end of a Channel
               
                
                % save for results this whole movie           
%                 save([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'],...
%                      'NA_feature_thisMD',...
%                     'Identifier_thisMD',...
%                     'CFMP_feature_ordered_MD_cell','CFMP_feature_MD_cell',...
%                     'ChMP_feature_MD_cell');

                  save([MD_ROOT_DIR,filesep,'movieData_NA_results_gathered.mat'],...
                    'Identifier_thisMD',...
                    'CFMP_feature_ordered_thisMD',...
                    'CFMP_feature_thisMD',...
                    'ChMP_feature_thisMD');

            end % end of if previous gathering exists for MD
            
            Identifier_thisML{1,iMD} = Identifier_thisMD;
            CFMP_feature_thisML{1,iMD} = CFMP_feature_thisMD;
            ChMP_feature_thisML{1,iMD} = ChMP_feature_thisMD;
            CFMP_feature_ordered_thisML{1,iMD} =CFMP_feature_ordered_thisMD;        
            
        end  % end of a MD        
        
         save([MD_ROOT_DIR,filesep,'movieList_NA_results_gathered.mat'],...
                    'Identifier_thisML',...
                    'CFMP_feature_ordered_thisML',...
                    'CFMP_feature_thisML',...
                    'ChMP_feature_thisML');
    end     % end of if previous gathering exists FOR a ML    
    
    
    %% network analysis plotting for vim screen
    plate_network_feature_plotting(Identifier_thisML,...
                    CFMP_feature_ordered_thisML,...
                    CFMP_feature_thisML,...
                    ChMP_feature_thisML);
    
    
    
end % end of a ML

