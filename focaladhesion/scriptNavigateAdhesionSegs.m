%% initialization
% load('/mnt/nas/CollagenV/ColV_Traction/Traction/092325-TFMAnalysis/ColI_2Kpa/movieList.mat')
load('/mnt/nas/CollagenV/ColV_Traction/Traction/092325-TFMAnalysis/ColV-2kPa/movieList.mat')
i=0;
%%  looping through
close
i = i+1;
curPath = fileparts(ML.movieDataFile_{i});
openfig([curPath filesep 'Adhesion Quantification' filesep 'imgs' filesep ...
    'figs' filesep 'imgNAFA001.fig'])

%% output
disp(curPath)