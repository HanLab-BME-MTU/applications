
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
MD.pixelSizeZ_ = 216;
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

%% Projection without ROI (full volume) in the lab frame of reference
disp('printMIP');
tic
printMIP(MDCrop);
toc;
disp('projectDynROI equivalent');
MDCrop.reset();
tic;
processProj1 = ProjectDynROIProcess(MDCrop,'Full MIP');
projectDynROI(MDCrop,'processSingleProj', processProj1);
toc; 
figure();
maxXY=processProj1.loadFrame(1,1);
imshow(maxXY);


% try using MD process
MDCrop.addProcess(processProj1);


%% Using a ROI in the lab frame of reference
% Building the ROI from tracks (Here the spindle poles)
rehash
[~, tracks] = detectPoles(MDCrop,'isoOutput',true);
ROI = tracks;
rehash
ref = FrameOfRef().genCanonicalRef(NumFrames);

tic;
processProj2 = ProjectDynROIProcess(MDCrop,'ROI-lab-ref');
% TODO: build a a projDynTOI Correctly inslide (setter/getter)
projectDynROI(MDCrop,ROI,ROI,'name', ['testDynProj'],'processSingleProj', processProj2, 'fringeWidth', 100,'insetFringeWidth',90);
toc;
figure();
maxXY=processProj2.loadFrame(1,1);
imshow(maxXY);

% attempt with adding process.
MDCrop.addProcess(processProj2);



%% Using a ROI and Ref on channel 2 only
rehash
ref = buildRefsFromTracks(tracks(1), tracks(2));
tic;
% Here "ROI-ref" specify the resulting output dir
processProj3 = ProjectDynROIProcess(MDCrop,'ROI-ref');
% In this object the output dir will be defined relatively to projessProj3
processRenderer = ProjectDynROIRendering(processProj3,'merged');

projectDynROI(MDCrop, ROI,'FoF',ref,'renderedChannel',2, ... 
    'intMinPrctil',[25 25],'intMaxPrctil',[99.5 99.99],...
    'processSingleProj',processProj3,'processRenderer',processRenderer, ...
    'fringeWidth',100,'insetFringeWidth',90,'processFrame',1:2);
toc;

% Create animation
anim=ProjAnimation(processRenderer,'ortho');
figure();imshow(anim.loadView(1));


%% Tracking KT and overlay detection and tracks (the same function are called in the trackKT function for debugging purposes)
%% The main issue of this implementation is that most of the required input 
%% should be implemented in a <DynROI> class. For example, per default, calling <overlayProjDetectionMovie>
%% on a set of detections should take a dynROI as input and set the detection in the proper FoF and 
%% mapping the object in the dynROI of interest. 

% Relaunch the ref detection so that "trackKT" can access the data.
buildAndProjectSpindleRef(MDCrop,'package',[]);
trackKT(MDCrop);
KTpack=MDCrop.searchPackageName('trackKT','selectIdx','last');

% load detection and register in Spindle FoF
tmp=load(KTpack.getProcess(1).outFilePaths_{2}); 
detection=tmp.movieInfo;
oDetections=Detections(detection);

%% The overlay require a processRenderer that actually process and merge channel
disp('Overlay KT detections');tic;
rehash;
processAllDetectOverlay = ProjectDynROIRendering(processRenderer,'allDetections');

myColormap=255*jet(256);
colorIndx=arrayfun(@(c) ceil(255*mat2gray(c.zCoord(:,1),[1,50]))+1,oDetections,'unif',0);
overlayProjDetectionMovie(processRenderer,'detections', oDetections , ... 
    'colorIndx',colorIndx, ...
    'colormap',myColormap,'name',['allDetections'],'process',processAllDetectOverlay);
figure();[~,~,~,ortho]=processAllDetectOverlay.loadFrame(1,1);imshow(ortho);
toc;

% load tracks and display
tmp=load(KTpack.getProcess(3).outFilePaths_{2}); 
kinTracksISO=TracksHandle(tmp.tracksFinal);
processTracksOverlay = ProjectDynROIRendering(processRenderer,'allTracks');

overlayProjTracksMovie(processAllDetectOverlay,'tracks', kinTracksISO, ... 
    'colorIndx',ceil(255*mat2gray([kinTracksISO.lifetime]',[1 MDCrop.nFrames_]))+1,'dragonTail',10, ...
    'colormap',myColormap,'name',['allTracks'],'process',processTracksOverlay);
figure();[~,~,~,ortho]=processTracksOverlay.loadFrame(1,2);imshow(ortho);
toc;



%% Test project1D for backward comp

% Build more complexe ROI and ref for fun

KT=kinTracksISO(1);
refPPKT=copy(ref);
refPPKT.genBaseFromZ(KT);
KTRefPPKT=refPPKT.applyBase(KT,'');
KTRefPPKTOpposite=KTRefPPKT.copy();
KTRefPPKTOpposite.x=-KTRefPPKTOpposite.x;

process_PPKSlice=ProjectDynROIProcess(MDCrop,'PPKTSlice');

project1D(  MD,[],'FoF',refPPKT,'dynPoligonREF',[refPPKT.applyBase([ROI],'') KTRefPPKT   KTRefPPKTOpposite], ...
    'name','PPKSlice','channelRender','greenRed','saveSingleProj',true, 'fringeWidth',[40,5,40], ...
    'processSingleProj',process_PPKSlice, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);

processAllDetectOverlay = ProjectDynROIRendering(process_PPKSlice,'allDetections');

myColormap=255*jet(256);
colorIndx=arrayfun(@(c) ceil(255*mat2gray(c.zCoord(:,1),[1,50]))+1,oDetections,'unif',0);
overlayProjDetectionMovie(process_PPKSlice,'detections', oDetections , ... 
    'colorIndx',colorIndx, ...
    'colormap',myColormap,'name',['allDetections'],'process',processAllDetectOverlay);
figure();[~,~,~,ortho]=processAllDetectOverlay.loadFrame(1,1);imshow(ortho);



%%
% Here "ROI-ref" specify the resulting output dir
processProj4 = ProjectDynROIProcess(MDCrop,'PPKSlice');
% In this object the output dir will be defined relatively to projessProj3
processRenderer = ProjectDynROIRendering(processProj4,'merged');

projectDynROI(MDCrop,[],refPPKT.applyBase([ROI KT],''),'FoF',refPPKT,'renderedChannel',2, ... 
    'intMinPrctil',[25 25],'intMaxPrctil',[99.5 99.99],...
    'processSingleProj',processProj4,'processRenderer',processRenderer, ...
    'fringeWidth',5,'insetFringeWidth',5,'processFrame',1:2);

processAllDetectOverlay = ProjectDynROIRendering(processRenderer,'allDetections');

myColormap=255*jet(256);
colorIndx=arrayfun(@(c) ceil(255*mat2gray(c.zCoord(:,1),[1,50]))+1,oDetections,'unif',0);
overlayProjDetectionMovie(processRenderer,'detections', oDetections , ... 
    'colorIndx',colorIndx,'radius',2, ...
    'colormap',myColormap,'name',['allDetections'],'process',processAllDetectOverlay);
figure();[~,~,~,ortho]=processAllDetectOverlay.loadFrame(1,1);imshow(ortho);
