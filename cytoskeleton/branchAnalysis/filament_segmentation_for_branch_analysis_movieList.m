function filament_segmentation_for_branch_analysis_movieList(ML, startML,Parameter_MD, varargin)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function

ip = inputParser;

ip.addRequired('ML', @(x) isa(x,'MovieList'));
ip.addRequired('startML',  @isnumeric);
ip.addRequired('Parameter_MD',@(x) isa(x,'MovieData') | isempty(x));

ip.addOptional('run_with_new_param',  0, @isnumeric);
ip.addOptional('input_parameter_set', [], @iscell);
ip.addOptional('save_old_data_tag',  [], @ischar);
ip.addOptional('whole_movie_filename', [], @ischar);


ip.parse(ML,startML,Parameter_MD,varargin{:});
whole_movie_filename= ip.Results.whole_movie_filename;
save_old_data_tag  =ip.Results.save_old_data_tag;
run_with_new_param = ip.Results.run_with_new_param;
input_parameter_set = ip.Results.input_parameter_set;

% the number of movies
movieNumber =  length(ML.movieDataFile_);

for iM  = startML:movieNumber
    
    clearvars -except 'movieNumber' 'iM' 'ML' 'whole_movie_filename' 'save_old_data_tag' 'run_with_new_param' 'input_parameter_set' 'Parameter_MD'
    
    close all;
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
    ticIDMD = tic;
    branchanalysis_load_MD_run_filament_analysis_package(MD, Parameter_MD, ...
        'save_old_data_tag', save_old_data_tag, 'run_with_new_param',run_with_new_param);
    display('Filament segmentation costed time for this movie:');
    toc(ticIDMD);
    
end