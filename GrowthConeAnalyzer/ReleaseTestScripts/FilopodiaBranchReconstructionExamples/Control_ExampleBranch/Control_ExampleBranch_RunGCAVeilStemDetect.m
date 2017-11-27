function Control_ExampleBranch_RunGCAVeilStemDetect( varargin )
% Control_ExampleBranch_RunGCAVeilStemDetect 
% Runs Steps I-IV for of GCA Segmentation (Veil/Stem Detection) for a 
% Branching Control Growth Cone 
% Images Include that used in Movie S3 (lower panel)
% Movie is truncated from original version for this software example to save
% time/space. 
%%
ip= inputParser;

ip.CaseSensitive = false;

ip.addParameter('OutputDirectory',[],@(x) ischar(x));
ip.addParameter('InputDirectory',[],@(x) ischar(x));
p = ip.Results;

ip.parse(varargin{:});

imageDirectory = ip.Results.InputDirectory;
outputDirectory = ip.Results.OutputDirectory;
if isempty(imageDirectory)
    imageDirectory =  ['/project/bioinformatics/Danuser_lab/shared/Maria/GCA/'...
        'ExampleData/FilopodiaBranchReconstructionExamples/Control_ExampleBranch/raw/Channels/Channel_1'] ;
end

if isempty(outputDirectory)
    outputDirectory =  ['/project/bioinformatics/Danuser_lab/shared/Maria/GCA/'...
        'ExampleData/FilopodiaBranchReconstructionExamples/Control_ExampleBranch/analysis'] ;
end
%% Make the movieData input/output object 
channel = Channel(imageDirectory);

% Constructor needs an array of channels and an output directory (for analysis)
MD = MovieData(channel,outputDirectory);

% Set the path where to store the MovieData object.
MD.setPath(outputDirectory);
MD.setFilename('Control_ExampleBranch_movieData.mat');

% Run sanityCheck on MovieData.
% Check image size and number of frames are consistent.
% Save the movie if successfull
MD.sanityCheck; % 

% Set some additional movie properties
MD.pixelSize_=216; % in nm after binning
MD.timeInterval_=5;% in sec
% Save the movieData
MD.save;

%% STEP I (Orange Block Figure S2)
% Run the first frame with the TSOverlays (Troubleshoot Overlays) for an
% example. 
GCAgetNeuriteOrientMovie(MD,'startFrame',1,'endFrame',1); 

% Run the rest of the movie with out the overlays 
GCAgetNeuriteOrientMovie(MD,'startFrame',2,'TSOverlays',false); 

% Notes for MD: 
% Create 'DetectStem' Process
% funName_ @GCAgetNeuriteOrientMovie
% funParams_ : p structure 
% inFilePaths_ : image channel 
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\I_neurite_orientation\Channel_1']
%  Saves: backboneInfo.mat (stemInfo for entire movie), params.mat, and some troubleshooting
%  visualizations if user desires. 
%% STEP II (Orange Block Figure S2)

GCAneuriteOrientConsistencyCheckMovie(MD) 
% Notes for MD: 
% Create 'CheckStem' Process
% funName_ @GCAneuriteOrientConsistencyCheckMovie.m
% funParams_ : p structure 
% inFilePaths_ : output file from Step I : loads backboneInfo.mat 
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\II_neurite_orientation_refinements\Channel_1']
%  Saves: updated backboneInfoFix.mat, params.mat, and some troubleshooting
%  visualizations if user desires/it is applicable. 
%% STEP III (Blue Block Figure S2) 
GCAReconstructVeilStemMovie(MD,'LocalThresholdPatchSize',40) 

% Notes for MD: 
% Create 'ReconstructVeilStem' Process
% funName_ @GCAReconstructVeilStemMovie.m
% funParams_ : p structure 
% inFilePaths_ : output file from Step II: loads backboneInfoFix.mat 
% (optional) can load in a mask file directory instead of the local thresholding
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\III_veilStem_reconstruction\Channel_1']
%  Saves: veilStem.mat, params.mat, and some troubleshooting
%  visualizations if user desires/it is applicable. 

%% STEP IV (Green Block Figure S2) 

GCAfindVeilStemLongestPathMovie(MD)  
% Notes for MD: 
% Create 'NeuriteLength' Process
% funName_ @GCAfindVeilStemLongestPathMovie.m
% funParams_ : p structure 
% inFilePaths_ : output file from Step III: loads veilStem.mat 
% outFilePaths_ : currently saved with p (parameter structure)
%  Default is [MD.outputDirectory_ filesep
%  'SegmentationPackage\StepsToReconstruct\IV_veilStem_length\Channel_1'
%  Saves: veilStem.mat (with additional fields, .neuriteLongPathIndices,.endPointLeadingProt),
%  neuriteLength.mat, params.mat, and some troubleshooting
%  visualizations if user desires/it is applicable. 

%% STEP V (OPTIONAL) Not applicable here. (Black Block Figure S2)

%% NOTE: After Steps I-IV(V) Veil/Stem Reconstruction complete: 
% Currently do a switch of the masks out of a global thresholding process 
% and then put into the protrusion software. 
%% Protrusion Vector Calc (For Veil/Stem Dynamics and Filopodia Orientation Calcs)



end

