function this_MD = vimscreen_load_ML_run_filament_analysis_package(ML, Parameter_MD, varargin)
% Function of single image filament segmentation with input this_MD from other
%               successfully segmented movie for the parameters in this_MD

% Input:      this_MD:              the movieData for data to be segmented
%             Parameter_MD:    a loaded this_MD with good segmentation parameters
%                              if none, put [], so a default setting will be used
%                              with (1) Otsu with smoothing 1 (2) mask
%                              refine with 1 object, (3) image flatten with
%                              square, (4) steerable filter with [1 2], (5)
%                              segmentation with geo based alpha=2
%             whole_movie_filename: optional input, if given, the whole
%                                movie statistics will be loaded from this filename, which
%                                will be used in the filament segmentsation; by default,
%                                nothing is input, so whole movie stat will be calculated,
%                                well, in this case, without any meaningful effect

% Created 2014.12 by Liya Ding, Matlab R2012b

% the number of movies


ip = inputParser;

ip.addRequired('ML', @(x) isa(x,'MovieList'));
ip.addRequired('Parameter_MD',@(x) isa(x,'MovieData') | isempty(x));

ip.addOptional('run_with_new_param',  0, @isnumeric);
ip.addOptional('input_parameter_set', [], @(x) iscell(x) | isempty(x));
ip.addOptional('save_old_data_tag',  [], @(x) ischar(x) | isempty(x));
ip.addOptional('whole_movie_filename', [], @(x) ischar(x) | isempty(x));

ip.parse(ML,Parameter_MD,varargin{:});

whole_movie_filename= ip.Results.whole_movie_filename;
save_old_data_tag  =ip.Results.save_old_data_tag;
run_with_new_param = ip.Results.run_with_new_param;
input_parameter_set = ip.Results.input_parameter_set;

movieNumber =  length(ML.movieDataFile_);

for iM  = 1:movieNumber
    
    clearvars -except 'movieNumber' 'iM' 'ML' 'Parameter_MD' 'varargin' 'run_with_new_param' ...
        'input_parameter_set' 'save_old_data_tag' 'whole_movie_filename'
        
    close all;
    
    % load this movie
    MD = MovieData.load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
                     
    vimscreen_load_MD_run_filament_analysis_package(MD,  Parameter_MD,...
         'run_with_new_param',run_with_new_param, ...
        'input_parameter_set',input_parameter_set,...
        'save_old_data_tag',save_old_data_tag,...
        'whole_movie_filename',whole_movie_filename);
end
