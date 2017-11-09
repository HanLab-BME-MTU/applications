function BA_output_GMST_ML_cell = branch_analysis_GMST_movieList(ML,half_size,min_branch_size_Threshold,filament_stat_flag,figure_flag,min_filament_length)
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

if(nargin<6)
    min_filament_length=0;
end

% the number of movies
movieNumber =  length(ML.movieDataFile_);
BA_output_GMST_ML_cell= cell(1,1);

for iM  = 1:movieNumber
    
    clearvars -except 'movieNumber' 'BA_output_GMST_ML_cell' 'iM' 'ML' 'figure_flag' 'half_size' 'min_branch_size_Threshold' 'filament_stat_flag' 'BA_output_GMST_ML_cell' 'min_filament_length'
    
    close all;
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
    tic   
    try
         BA_output_GMST_ML_cell{1, iM} = branch_analysis_GMST_movieData(MD, half_size, min_branch_size_Threshold, filament_stat_flag, figure_flag,min_filament_length); 
    catch
         display(['iM:', num2str(iM),', corrupted.']);
        continue;
    end
    display(['iM:', num2str(iM),', costed:']);
    toc
end

 ML_ROOT_DIR = ML.outputDirectory_;
 save([ML_ROOT_DIR,filesep,'movieList_BA_output_GMSTminLength.mat'],'BA_output_GMST_ML_cell');

 
 try
     plot(A);
 end