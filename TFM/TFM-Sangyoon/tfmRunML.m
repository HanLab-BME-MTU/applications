function [] = tfmRunML(MLPath,useParFor)
if nargin<2
    useParFor = false;
end
%% set up
% analysisFolder = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2017-06-29//ChoK1_shRNA_WT_Rescue_FACS_5kPa_006';
if isa(MLPath,'MovieList')
    ML=MovieList.load(MLPath.getFullPath); %,'askUser',false,'askUserChannel',false);
else
    ML=MovieList.load(MLPath,'askUser',false); %,'askUserChannel',false);
end

nMovies = numel(ML.movieDataFile_);
MDAll = ML.movies_;
if useParFor
    parfor ii=1:nMovies
        curMD = MDAll{ii};
        tfmRun(curMD)
    end
else
    for ii=1:nMovies
        curMD = MDAll{ii};
        tfmRun(curMD)
    end
end
ML.save
end