function save_MIP_features_directions_coding(mov_location,filename,outname,params_location,add_path)

%% ================================================================================
% save_MIP_features_directions_coding
% -------------------------------------------------
% Runs MIP coding for all directions on one movie.
%
% Inputs:
% ~~~~~~
% mov_location - path to video (usually: InputDir)
% filename - specific video file name
% outname (full path) - out file name for the descriptor  (suffix .mat!)
% params_location (full path) - path to the currecnt test params 
% add_path - path of the code. from GUI dir:  .. (i.e. MIPCode/Code)
%
%
% Copyright 2012, Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
%
%% ================================================================================

% Add code path:
if(~iscell(add_path))
    addpath(genpath(add_path));
else
    for i=1:length(add_path)
        addpath(genpath(add_path{i}));
    end
end

% Load params:
load(params_location);
MIP_params = params.MIP_params;

% Get the directions for coding
directions = MIP_params.MIP_code.batch_run_directions;

% Define the features' root path 
org_outname = outname;
[out_path,mov_name,ext] = fileparts(org_outname);
ind = strfind(out_path,params.features.dir_name);
features_parent_path = out_path(1:ind-1);

% Alignemnt:
% --------------

% before going over directions (for computational reasons)
% if the direct falg is on - save transformations and valid mask 
% on the features root directory - reload them when needed.
if(MIP_params.align.en && MIP_params.align.precompute.en)
    if(isempty(MIP_params.align.precompute.dir_name))
         % Define the dir for precomputed wrap - relative to features path
        MIP_params.align.precompute.dir_name = fullfile(features_parent_path,'wrapped_frames');
    end
    % define out wrap frame of specific file
    precomp_full_name = fullfile(MIP_params.align.precompute.dir_name ,[filename(1:end-4),'_wrapped.mat']);    
    if(~exist(precomp_full_name,'file'))
        % Debug dir location - if empty:
        if(MIP_params.align.dbg_out && (isfield(MIP_params.align,'dbg_dir') || isempty(MIP_params.align.dbg_dir)))
            MIP_params.align.dbg_dir = fullfile(MIP_params.align.precompute.dir_name,'debug_data');
        end
        % call for computing alignemnt
        MIP_params.MIP_code.code_shift = 0;
        MIP_params.MIP_dbg_dumps = false;%dont out debug encode for precomputing alignmrnt.
        [P,F] = get_MIP_points(mov_location,filename,MIP_params);
        MIP_params.MIP_code.code_shift  = params.MIP_params.MIP_code.code_shift; 
        MIP_params.MIP_dbg_dumps = params.MIP_params.MIP_dbg_dumps;
    end
end
MIP_params.align.dbg_out = false;%assume here that precomputed alignment.

% Run MIP and save to file:
% --------------------------------

% run directions , with changing code_shift:
if(~isempty(directions))   
    
    for d = 1:length(directions)           
        % get per direction outname:
        cur_out_path = fullfile(features_parent_path,['direction_',num2str(directions(d))]);
        outname = fullfile(cur_out_path,params.features.dir_name,[mov_name,ext]);
        if(~exist(fullfile(cur_out_path,params.features.dir_name),'dir'))
            mkdir(fullfile(cur_out_path,params.features.dir_name));
        end
        % if exist dont code again, unless forced to overide
        if (exist(outname,'file') && (~params.features.overide))
            continue;
        end
        if(d==1)% for debug only...
            fprintf(1,'Extracting Features mov : %s\n',filename);
        end
        % Debug dir location:('params...' doesnt change during iterations)
        if(MIP_params.MIP_dbg_dumps && (~isfield(params.MIP_params,'dbg_dir') || isempty(params.MIP_params.dbg_dir)))
            MIP_params.dbg_dir = fullfile(cur_out_path,'Debug_data');
        end
        
        % The actual coding:
        MIP_params.MIP_code.code_shift = directions(d);
        [P,F] = get_MIP_points(mov_location,filename,MIP_params);
        F = sparse(double(F));
        P = single(P);
    
        if(isempty(F))
            warning(['Movie ',outname,' has an empty coding results']);
        end
        % Save result to file:
        save(outname,'F','P','-v7.3');
        fprintf(1,'Mov : %s : Direction %i Done. \n',filename,directions(d));
    end
    
else
    % No multiple directions, not supported
    error('No directions set, Please set the vector ''batch_run_directions''');
end

end % of function: save_MIP_features_directions_coding


