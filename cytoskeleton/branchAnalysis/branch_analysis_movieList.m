function BA_output_ML_cell = branch_analysis_movieList(ML,half_size,min_branch_size_Threshold,filament_stat_flag,figure_flag)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function


if(nargin<2)
    half_size=150;
end

if(nargin<3)
    min_branch_size_Threshold=100;
end

if(nargin<4)
    filament_stat_flag=0;
end

if(nargin<5)
    figure_flag=0;
end

% the number of movies
movieNumber =  length(ML.movieDataFile_);
BA_output_ML_cell= cell(1,1);

for iM  = 3:movieNumber
    
    clearvars -except 'movieNumber' 'BA_output_ML_cell' 'iM' 'ML' 'figure_flag' 'half_size' 'min_branch_size_Threshold' 'filament_stat_flag'
    
    close all;
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
                     
    BA_output_ML_cell{1, iM} = branch_analysis_movieData(MD, half_size, min_branch_size_Threshold, filament_stat_flag, figure_flag);  
    
end

 ML_ROOT_DIR = ML.outputDirectory_;
 save([ML_ROOT_DIR,filesep,'movieList_BA_output_balloon.mat'],'BA_output_ML_cell');

 branch_analysis_moiveList_results_gather(ML);
 