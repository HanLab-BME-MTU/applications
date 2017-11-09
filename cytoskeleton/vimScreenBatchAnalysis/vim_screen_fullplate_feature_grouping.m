function vim_screen_fullplate_feature_grouping(plate_analysis_folder,feature_index, results_saving_folder)
% Function for gathering and organizing the network analysis resulting features for the full plate of 384 wells.

% Liya Ding
% Matlab 2015a
% 2015.12

% Note: due to the size of the data and analysis, the raw images and some
% results will be removed after the FilamentAnalysisPackage and network
% analysis job finishes. So the movieData structure is broken. To deal with
% this, given the vim screen struture is fixed, the folders here are hard coded.
% For the future, if the structure updates, the code need to be updated as well.

% Input: plate_analysis_folder :              The root folder of the plate analysis. (String of the folder name).
%        feature_index:                       The features interested. If not given, use all 28 possible features.
%                                               ( note that after re-organization, there are 33 features.)
%        results_saving_folder:               The folder to save all the results for gathered features.
%                                               If not given, this will be the movieList folder.



% if no feature index is given, use all features
if(nargin<2)
    feature_index=ones(28,1);
end

% if no output saving folder is given, use the movieList folder under plate analysis folder

if(nargin<10)
    results_saving_folder = [plate_analysis_folder, filesep,'row_col_movieList'];
end

movieNumber = 384; % fixed

ML_name = [plate_analysis_folder, filesep,'row_col_movieList',filesep,'movieList.mat'];

% load directly, not movieList.load since the list of MD are broken.
load(ML_name); % then we have ML

movieNumber = numel(ML.movieDataFile_);

% we have two channels, set cell of two for the whole ML
NA_feature_thisML = cell(2,1);
Identifier_thisML = cell(2,1);
NA_feature_whole_thisML = cell(2,1);

% for each well(MD)
for iMD  = 1 : movieNumber
    
    % load one by one
    NA_feature_thisMD = cell(2,1);
    Identifier_thisMD = cell(2,1);
    NA_feature_whole_thisMD = cell(2,1);
    
    
    try
        % load MD directly, not movieData.load since the MD is broken.
        load(ML.movieDataFile_{iMD});
    catch
        disp('The MD file is missing');
        continue;
    end
    
    % Now MD in workspace
    % check index
    display_msg_flag=0;
    package_process_ind_script;
    
    nChannel = numel(MD.channels_);
    
    nFrame = MD.nFrames_;
    
    % check channel 1 and 2, for vim and mt, fixed
    for iChannel = 1 : 2
        display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);
        
        % hard coded folders
        MD_channel_analysis_folder = [MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel',num2str(iChannel),filesep,'analysis_results'];
        
        % try to given the image names
        % when the raw images are removed, the MD doesn't have the correct image folders
        % then use the movie index directly to get the resulting analysis
        % filenames, assuming the order of the wells are correct
        try
            Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
            filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
        catch
            % assuming 16 rows, 24 cols
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
        if(exist(MD_channel_analysis_folder,'dir'))
            % if it exist, load the existing analysis
            
            for iFrame = 1 : nFrame
                display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
                load_ch_frame_flag =0;
                
                % the two possible mat file names due to code verions with
                % new and old code discrepency
                
                if(exist([MD_channel_analysis_folder,filesep,'network_analysis_feature_ch_',num2str(iChannel),...
                        '_frame_',num2str(iFrame),'.mat'],'file'))
                    load([MD_channel_analysis_folder,filesep,'network_analysis_feature_ch_',num2str(iChannel),...
                        '_frame_',num2str(iFrame),'.mat'],'output_feature');
                    load_ch_frame_flag=1;
                end
                
                if(exist([MD_channel_analysis_folder,filesep,'network_feature_ch_',num2str(iChannel),'_',...
                        filename_short_strs{iFrame},'.mat'],'file'))
                    load([MD_channel_analysis_folder,filesep,'network_feature_ch_',num2str(iChannel),'_',...
                        filename_short_strs{iFrame},'.mat'], 'output_feature');
                    load_ch_frame_flag = 1;
                end
                
                if(load_ch_frame_flag==0)
                    continue;
                end
                
                % arrange the features, all into column vectors(density, keep only ROI part)
                % also organize the
                output_feature_arranged = arrange_network_features(output_feature);
                
                % organize profiles for the P1 P2 ratio analysis
                % here the features change from 28 to 33
                % look into the function organize_profiles to check more of
                % the feature index information
                
                output_feature_reorganized = organize_profiles(output_feature_arranged,feature_index);
                
                % update the index as well
                feature_index_reorganized = feature_index;
                
                % the original feature index 19 is for the profiles, update accordingly
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
                
                % short script to put all the sites of this movie into the the movie pools
                % Not using functions since the variables are big sized
                put_frame_pool_into_MD_pool; % after this we have NA_feature_whole_thisMD
                
                
            end% end of a frame
            
        end % end of if the analysis for output folder exist
        
        % since each movie is each well, get the well percentiles, (for a channel)
        ChMP_feature_thisMD{iChannel} = extract_mean_percentiles_well(NA_feature_whole_thisMD, iChannel, valid_requested_feature_index);
        
    end % end of a Channel
    
    
    save([results_saving_folder,filesep,Identifier_thisMD{iChannel, iFrame},'_movieData_well_NA_results_pool_gathered.mat'],...
        'Identifier_thisMD',...
        'ChMP_feature_thisMD',...
        'NA_feature_whole_thisMD',...
        'valid_feature_index',...
        'valid_requested_feature_index');
    
    Identifier_thisML{1,iMD} = Identifier_thisMD;
    ChMP_feature_thisML{1,iMD} = ChMP_feature_thisMD;
    
end  % end of a MD

save([results_saving_folder,filesep,'movieList_plate_NA_results_gathered.mat'],...
    'Identifier_thisML',...
    'ChMP_feature_thisML');
