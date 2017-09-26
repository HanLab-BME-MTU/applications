
%% Building test dataset

% MATLAB-based testing/performance suite for Utrack3D
% Andrew R. Jamieson 2017
% Test Utrack3D+tips Package

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Preconditions
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_paths = path;
start_dir = pwd;


disp(['PWD:', pwd]);
% Dump path for debugging
s_path = strsplit(path,':');
s_match = (cellfun(@(x) regexp(x,'toolbox'), s_path, 'UniformOutput', false))';
matlab_paths = s_path(cellfun(@isempty, s_match))';
disp('    [MATLAB] current top-level paths....');
disp(matlab_paths);

disp(['Java heap max: ' num2str(java.lang.Runtime.getRuntime.maxMemory/1e9) 'GB'])
disp('Starting Utrack3D package test script -- 2 Channel EB1 - HeLa');

%----Initialization of temp dir
package_name = 'UTrackPackage3D';
t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
mkdir(tmpdir);
cd(tmpdir);


% Download test data for u-track (Use BioHPC internal cloud account danuserweb)
url = 'https://lamella.biohpc.swmed.edu/index.php/s/RbW7OPNuPiyEndK/download'
zipPath = fullfile(tmpdir, 'EB1_3D.zip');

try 
    urlwrite(url, zipPath);
catch
    disp('lamella failed, try using local file cache..')
    copyfile('C:\Users\Andrew\Data\raw\3D\Tracking\EB1_3D.zip', zipPath)
end

unzip(zipPath, tmpdir);


% Initialize MovieData ch1
ch0 = fullfile(tmpdir, 'ch0'); % Comets (use watershedApple)
ch1 = fullfile(tmpdir, 'ch1'); % kinetichores (use pointSrc3D)

% Analysis Output Directory
saveFolder = [tmpdir filesep 'Analysis'];
analysis_dir = saveFolder;
mkdir(saveFolder);


% Imaging/Microscope Parameters
NA = 1.42;
pixelSize = 80;
timeInterval = 1;

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Channels
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve current location
% Channel creation
% Create a channels object
Chan0 = Channel(ch0);
Chan1 = Channel(ch1);
channels = [Chan0 Chan1];

% channel.fluorophore_='alexa647';
% channel.emissionWavelength_=name2wavelength('alexa647')*1e9;
% channel.imageType_='TIRF';

% channel(2) = Channel(ch1);
% channel(2).fluorophore_='';
% channel(2).emissionWavelength_=name2wavelength('egfp')*1e9;
% channel(2).imageType_='';
%----------------------------------------

% MovieData creation
outputDir = analysis_dir; %fileparts(analysis_dir);
MD = MovieData(channels, outputDir);

movieDataFileName = 'movieDataScriptTestHeLaEB1_3D.mat';
MD.setPath(outputDir);
MD.setFilename(movieDataFileName);

% Imaging/Microscope/Movie Parameters
% MD.numAperture_= 1.49; 
MD.pixelSize_= 100;
MD.pixelSizeZ_ = 400;
MD.camBitdepth_= 16;
MD.timeInterval_ = 1;
MD.notes_= 'HeLa A1 EB1 tracking 3D testing movie with new Utrack3Dpackage CI - 2 channels'; 
MD.sanityCheck;
MD.save;


% Load the movie
clear MD;
MD = MovieData.load(fullfile(outputDir, movieDataFileName));
MD.reset();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% The attached script was located at Z:\repo\utsw-ssh\applications\mtTracker\EB3-3D\graphics\testScript
% The idea is that the class ProjectDynROIProcess inherits from ComputeMIPProcess , 
% when launching the movieViewer GUI we want to access the different processes. 
% Switching between the two, or even better, playing both projection at the same time.  
% In my hand, I couldn’t see both in the same context.
% Then we need to manage the detection and tracks in those associated dynamical ROI.
% Then we need to allow the specification of ROI from the GUI. We need to spec that after some more progress.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%testing on prometaphase data
% allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
% allMovieToAnalyse=allMovieToAnalyse(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell),:);
% outputPath='/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/proudot/u-track-3D/dynROI/testProjectSingle';
% 
% MD=MovieData.loadMatFile(allMovieToAnalyse.analPath{1});
NumFrames = 5;
MDCrop = crop3D(MD,(MD.getChannel(1).loadStack(1)),'keepFrame',1:NumFrames,'name','testDynProj');

%%%% <<< WRAP crop3D as process. ?

%% Basic MIP
disp('printMIP');
tic
printMIP(MDCrop);
toc;
disp('projectDynROI equivalent');
MDCrop.reset();
tic;
processProj1 = ProjectDynROIProcess(MDCrop);
projectDynROI(MDCrop,'processSingleProj', processProj1);
toc; 
figure();
imshow(imread(sprintfPath(processProj1.outFilePaths_{7},1)));

% try using MD process
MDCrop.addProcess(processProj1);


%% Using a ROI
processProj2 = ProjectDynROIProcess(MDCrop);
ref = FrameOfRef().genCanonicalRef(NumFrames);
[poleMovieInfo, tracks] = detectPoles(MDCrop,'isoOutput',true);
ROI = tracks;

processVolMask = ProjectDynROIProcess(MDCrop);

tic;
projectDynROI(MDCrop,ROI, ...
    'name', ['testDynProj'], ...
    'channelRender','grayRed','processSingleProj', processProj2, ...
    'processMaskVolume', processVolMask,'crop','manifold', ... 
    'intMinPrctil',[20 98], 'intMaxPrctil', [100 100], 'fringeWidth', 50,'insetFringeWidth',10);
toc;
figure();
imshow(imread(sprintfPath(processProj2.outFilePaths_{7},1)));

% attempt with adding process.
MDCrop.addProcess(processProj2);



%% Using a ROI and Ref
ref = buildRefsFromTracks(tracks(1), tracks(2));
tic;
processProj3 = ProjectDynROIProcess(MDCrop);
projectDynROI(MDCrop, ROI,'FoF',ref, ...
    'name',['testDynProj-no-mask'],...%'suppressROIBorder',true, ...
    'channelRender','grayRed','processSingleProj',processProj3, ...
    'processMaskVolume',processVolMask,'crop','manifold', ... 
    'intMinPrctil',[20 98],'intMaxPrctil',[100 100],'fringeWidth',50,'insetFringeWidth',10);
toc;

figure();
imshow(imread(sprintfPath(processProj3.outFilePaths_{7},1)));
MDCrop.addProcess(processProj3);