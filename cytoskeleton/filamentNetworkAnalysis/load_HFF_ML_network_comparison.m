function load_HFF_ML_network_comparison(ML,iChannel1, start_frame1,...
    iChannel2, start_frame2, ...
    radius, show_save_everything_flag, ...
    longest_radius, sigma_gaussian, sigma_d, sigma_theta)
% function to compare two networks for a movieList

% function to do network comparison for a whole movielist
% Liya Ding, April, 2015
%
% Input:
%   ML:     The movieList object loaded before running this function
%           radius: the definition of neighborhood
%           figure_flag: 1 to plot histgrams, 0 not to
%           save_everything_flag: 1 to save the plot, 0 not to
%          iChannel1, iChannel2:
%                       The channel in the movie to be compared.
%           start_frame1, start_frame2:
%                       the start frame of comparison; if 1 and 1, then
%                       same time; if 1 and 2 then the dynamics
%          radius:      the local neighborhood size for network comparison
%          save_everything_flag:
%                       whether to save all the figure during detailed process
%                       usually set to 0, to keep only the useful ones

if(~exist('radius','var'))
    radius = 20;
end

if(~exist('longest_radius','var'))
    longest_radius = 2*radius;
end

if(~exist('figure_flag','var'))
    figure_flag = 0;
end

if(~exist('show_save_everything_flag','var'))
    show_save_everything_flag = 0;
end


flag_default=0;


if(~exist('sigma_gaussian','var') && ~exist('sigma_d','var') && ~exist('sigma_theta','var'))
   flag_default = 1;
end

if(~exist('sigma_gaussian','var'))
  sigma_gaussian = 3*radius/8;
    flag_default = 1;
end

if(~exist('sigma_d','var'))
  sigma_d = sqrt(3)*radius/4;
    flag_default = 1;
end


if(~exist('sigma_theta','var'))
 sigma_theta = pi/(2*sqrt(3));
    flag_default = 1;
end


%%
movieNumber =  length(ML.movieDataFile_);
similarity_scoremap_ML_cell = cell(1,movieNumber);
difference_map_ML_cell = cell(1,movieNumber);

for iM  = 1 :movieNumber
    
    clearvars -except 'movieNumber' ...
        'iM' 'ML' 'figure_flag' 'radius'...
        'show_save_everything_flag' ...
        'iChannel1' 'start_frame1' ...
        'iChannel2' 'start_frame2' ...
        'longest_radius' 'sigma_gaussian' 'sigma_d' 'sigma_theta'...
        'similarity_scoremap_ML_cell' 'difference_map_ML_cell'

    close all;
    
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
    
    [similarity_scoremap_cell,difference_map_cell] = load_HFF_MD_network_for_dynamics_compare(ML.movieDataFile_{iM},iChannel1, start_frame1,...
    iChannel2, start_frame2, ...
    radius,show_save_everything_flag,...
    longest_radius,sigma_gaussian, sigma_d, sigma_theta);

    similarity_scoremap_ML_cell{1,iM} = similarity_scoremap_cell;
    difference_map_ML_cell{1,iM} = difference_map_cell;
end

save([ML.outputDirectory_,filesep,'similarity_scoremap_ML_result_all.mat'],'similarity_scoremap_ML_cell','difference_map_ML_cell');

