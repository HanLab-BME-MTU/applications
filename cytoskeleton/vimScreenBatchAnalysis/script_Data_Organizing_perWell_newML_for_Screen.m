% Data organizing script for screen data for corruptted ML after
% relocation
 
Analysis_Dir = '/project/bioinformatics/Danuser_lab/vimscreen/analysis/lding/fromTony/P042_process_20151115';

%% Load the movieData object and save movieList for a  plate.

Movie_number = 16*24;

movie_cell = cell(1,Movie_number);

for iMovie = 1 : Movie_number
    
    col_num = mod(iMovie,24);
    
    row_num = floor(iMovie/24);
    
    if(col_num ==0)
        row_num = row_num-1; 
        col_num = 24;
    end
        
    row_col_ID = [char(row_num+65) '_' num2str(col_num,'%02d')];        
   
    % where to save it
    outputDirectory = [Analysis_Dir,filesep,'well_',row_col_ID];
     %% load MD
    MD = MovieData.load([outputDirectory,filesep,'movieData.mat']);   
    MD.sanityCheck();
    
    %% save movieData to movieList
    
    movie_cell{1,iMovie}= MD;
end

%% save movieList for the whole plate

mkdir([Analysis_Dir,filesep,'row_col_movieList']);

ML = MovieList(movie_cell,[Analysis_Dir,filesep,'row_col_movieList']);
ML.setPath([Analysis_Dir,filesep,'row_col_movieList']);
ML.setFilename('movieList.mat')
ML.sanityCheck();
