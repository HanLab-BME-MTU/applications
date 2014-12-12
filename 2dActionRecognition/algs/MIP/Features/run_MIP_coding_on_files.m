function features_dir = run_MIP_coding_on_files(test_dir,params,database_info,add_path)

%% ================================================================================
% run_MIP_coding_on_files
% ---------------------------------
% Run MIP coding for all requested files
% 
% Inputs:
% ~~~~~~
% test_dir - specific dir for the current test
% params - params structure created according to get_feature_params_runtime.m
% database_info - database_info is created based on the video files in the InputDir
% assumes the following structure of database_info: 
% --- database_info.location- start path of all movies (InputDir)
% --- database_info.db_data - {filename,relative path to location}
% add_path - path from GUI dir:  .. (i.e. MIPCode/Code)
%
% Outputs:
% ~~~~~~~
% features_dir -  where to write the output descriptors 
% (for future use - changes inside function)
%
% Copyright 2012, Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
%
%% ================================================================================

% Create dir for the current test
if(~exist(test_dir,'dir'))
    mkdir(test_dir);
end

% Change feature location to a given location - not neccesarily connected
% to the test (for example default output dir) 
if( isfield( params.features, 'location') && ~isempty( params.features.location) ) 
     features_dir = fullfile(params.features.location);
else
    % Define the features_dir relative to the current test_dir
    features_dir = fullfile(test_dir);
end
if(~exist(features_dir,'dir'))
    mkdir(features_dir);
end

% Next is for parallel use on windows:
if(params.flow.parfor_en)
    isOpen = matlabpool('size') > 0;
    if(~isOpen)
        matlabpool open;
    end
end
	
% Save params and database_info to the current run test dir:
params_location = fullfile(test_dir,'feature_params.mat');
save(params_location,'params');
save(fullfile(test_dir,'database_info.mat'),'database_info');

% InputDir
mov_location = database_info.location;


% Features extraction - loop on files:
% -------------------------------------------

% for parallel use on windows:
if(params.flow.parfor_en)		
    parfor m=1:size(database_info.db_data,1)
        filename = database_info.db_data{m,1};
        % fprintf(1,'%d, Excracting Features mov : %s\n',m,filename);
        outname = fullfile(features_dir,params.features.dir_name,[database_info.db_data{m,1}(1:end-4),'.mat']);
        cur_mov_path = fullfile(mov_location,database_info.db_data{m,2});
        save_MIP_features_directions_coding(cur_mov_path,filename,outname,params_location,add_path);
    end
% regular run / Unix
else
    for m=1:size(database_info.db_data,1)
        filename = database_info.db_data{m,1};
        % fprintf(1,'%d, Excracting Features mov : %s\n',m,filename);
        outname = fullfile(features_dir,params.features.dir_name,[database_info.db_data{m,1}(1:end-4),'.mat']);
        cur_mov_path = fullfile(mov_location,database_info.db_data{m,2});
        save_MIP_features_directions_coding(cur_mov_path,filename,outname,params_location,add_path);
    end		
end

end % of function: run_MIP_coding_on_files