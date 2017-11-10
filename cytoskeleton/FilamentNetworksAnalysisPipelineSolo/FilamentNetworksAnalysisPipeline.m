function FilamentNetworksAnalysisPipeline(vim_image_folder,mt_image_folder, analysis_folder)
% function for the whole vimentin-mt filament network interaction analysisprocess

% locate and build MD
MD = locate_two_channels_build_MD(vim_image_folder,mt_image_folder, analysis_folder);

% filament network segmentation
MD = load_MD_run_filament_segmentation(MD,[]);

% network feature pooling
load_MD_network_for_analysis(MD);

% similarity analysis
load_2_MD_network_for_dynamics_compare([analysis_folder,filesep,'movieData.mat'],1,1,...
[analysis_folder,filesep,'movieData.mat'],2,1,20,1);

% dynamic(self similarity) analysis
load_2_MD_network_for_dynamics_compare([analysis_folder,filesep,'movieData.mat'], 1,1,...
[analysis_folder,filesep,'movieData.mat'],1,2,20,1);

load_2_MD_network_for_dynamics_compare([analysis_folder,filesep,'movieData.mat'], 2,1,...
[analysis_folder,filesep,'movieData.mat'],2,2,20,1);
