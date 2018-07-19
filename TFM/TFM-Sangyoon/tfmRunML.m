function [] = tfmRunML(ML)
%% set up
% analysisFolder = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2017-06-29//ChoK1_shRNA_WT_Rescue_FACS_5kPa_006';
ML=MovieList.load([ML.movieListPath_ filesep ML.movieListFileName_]);
nMovies = numel(ML.movies_);
for ii=1:nMovies
    curMD = ML.movies_{ii};
    tfmRun([curMD.movieDataPath_ filesep curMD.movieDataFileName_])
end
ML.save
