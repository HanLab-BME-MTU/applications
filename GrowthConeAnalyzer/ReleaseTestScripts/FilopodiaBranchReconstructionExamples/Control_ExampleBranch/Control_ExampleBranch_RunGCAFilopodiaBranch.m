function [ output_args ] = Control_ExampleBranch_RunGCAFilopodiaBranch( varargin )
% Control_ExampleBranch_RunGCAVeilStemDetection 

ip= inputParser;

ip.CaseSensitive = false;

ip.addParameter('OutputDirectory',[],@(x) ischar(x));
%ip.addParameter('InputDirectory',[],@(x) ischar(x)); 
p = ip.Results;
 
ip.parse(varargin{:});
%imageDirectory = ip.Results.InputDirectory; 
outputDirectory = ip.Results.OutputDirectory;
% if isempty(imageDirectory) 
%     imageDirectory =  ['/project/bioinformatics/Danuser_lab/shared/Maria/GCA/'...
%         'ExampleData/FilopodiaBranchReconstructionExamples/Control_ExampleBranch/raw/Channels/Channel_1'] ; 
% end 

if isempty(outputDirectory)
    outputDirectory =  ['/project/bioinformatics/Danuser_lab/shared/Maria/GCA/'...
        'ExampleData/FilopodiaBranchReconstructionExamples/Control_ExampleBranch/analysis'] ; 
end 
    
 load([outputDirectory filesep 'Control_ExampleBranch_movieData.mat'])
%% STEP VI: (Purple Block in Figure S2) 
GCAReconstructFilopodiaMovie(MD); 
% Notes for MD: 
% Create 'FiloBranchReconstruct' Process
% funName_ @GCAReconstructFilopodiaMovie
% funParams_ : p structure 
% inFilePaths_ : output of Step IV (V-if run) of the GCA package 
% loads the veilStem.mat
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\VI_filopodiaBranch_reconstruction']
%  Saves: filoBranch.mat (filopodia and branch information for entire movie), params.mat, and some troubleshooting
%  visualizations if user desires. 
%% STEP VII (Pink Block in Figure S2)
% only make the troubleshoot overlays for the first frame 
GCAfitFilopodiaMovie(MD,'startFrame',1,'endFrame',1);

% run through the rest without making the overlays
GCAfitFilopodiaMovie(MD,'startFrame',2,'TSOverlays',false); 
% Notes for MD: 
% Create 'FilopodiaFit' Process
% funName_ @GCAfitFilopodiaMovie.m
% funParams_ : p structure 
% inFilePaths_ : output of Step VI of the GCA package 
% loads the filoBranch.mat
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\VII_filopodiaBranch_fits\Channel_1']
%  Saves: filoBranch.mat with updated fields (filopodia and branch information for entire movie), params.mat, and some troubleshooting
%  visualizations if user desires.
%%
%% Segmentation  Complete 
%%
%% Start Additional Functions for Measurements/Troubleshooting/Visualizations   
GCAAddFilopodiaNormalizedIntensityMovie(MD); % might add to above easier for user: depending on what filtering you require
% might not need it 
% Notes for MD: 
% Create 'AddFilopodiaIntensity' Process
% funName_ @GCAAddFilopodiaNormalizedIntensityMovie.m
% funParams_ : p structure 
% inFilePaths_ : output of Step VII of the GCA package 
% loads the filoBranch.mat
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\VII_filopodiaBranch_fits\Channel_1']
%  Saves: filoBranch.mat with updated fields (filopodia and branch information for entire movie), params.mat, and some troubleshooting
%  visualizations if user desires.
%% Troubleshoot FiloBranch Reconstruction Visualization 
% Make Video3_Fig1D (bottom panel) frame 7 
% 
% Make reconstruction movies for both frames 1 and 7 
GCATroubleshootMakeMovieOfReconstructMovie(MD,'frames',[1,7]); 
% 
% Not sure this really appropratiate for a process. It is more a visualization for 
% troubleshooting associated with a reconstruction Step VI. 
% inFilePaths_ : output of Step VI (for troubleshooting reconstruction before
% fitting) or Step VII (final- likewise adds the filopodia colored by
% length plots): 
% reads in filoBranch.mat
% Step IV from GCA (veil/stem detection): reads in veilStem.mat (could
% likewise been from StepIII,Step V).  
% outFilePaths_ :
%  Default is [MD.outputDirectory_ filesep
% Reconstruct_Movies]
