%% Building test dataset
%testing on prometaphase data
allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
allMovieToAnalyse=allMovieToAnalyse(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell),:);
outputPath='/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/proudot/u-track-3D/dynROI/testProjectSingle';

MD=MovieData.loadMatFile(allMovieToAnalyse.analPath{1});
MDCrop=crop3D(MD,(MD.getChannel(1).loadStack(1)),'keepFrame',1:5,'name','testDynProj');

% Basic MIP
disp('printMIP');
tic
printMIP(MDCrop);
toc;
disp('projectDynROI equivalent');
tic;
processProj=ProjectDynROIProcess(MDCrop);
projectDynROI(MDCrop,'processSingleProj',processProj);
toc; 
figure();
imshow(imread(sprintfPath(processProj.outFilePaths_{4},1)));

%% Using a ROI
processProj=ProjectDynROIProcess(MD);
ref=FrameOfRef().genCanonicalRef(5);
[poleMovieInfo,tracks]=detectPoles(MDCrop,'isoOutput',true);
ROI=tracks;

processVolMask=ProjectDynROIProcess(MD);
tic;
projectDynROI(MDCrop,ROI, ...
    'name',['testDynProj'], ...
    'channelRender','grayRed','processRenderer',processProj, ...
    'processMaskVolume',processVolMask,'crop','manifold', ... 
    'intMinPrctil',[20 98],'intMaxPrctil',[100 100],'fringeWidth',50,'insetFringeWidth',10);
toc;
figure();
imshow(imread(sprintfPath(processProj.outFilePaths_{4},1)));

%% Using a ROI and Ref
ref=buildRefsFromTracks(tracks(1),tracks(2));
tic;
projectDynROI(MDCrop,ROI, ...
    'name',['testDynProj-no-mask'],'suppressROIBorder',true, ...
    'channelRender','grayRed','processRenderer',processProj, ...
    'processMaskVolume',processVolMask,'crop','manifold', ... 
    'intMinPrctil',[20 98],'intMaxPrctil',[100 100],'fringeWidth',50,'insetFringeWidth',10);
toc;
figure();
imshow(imread(sprintfPath(processProj.outFilePaths_{4},1)));