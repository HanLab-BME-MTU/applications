%% Building test dataset
% testing on prometaphase data
allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
allMovieToAnalyse=allMovieToAnalyse(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell),:);
outputPath='/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/proudot/u-track-3D/dynROI/testProjectSpeed';

movieIndex=find(ismember(allMovieToAnalyse.Setup_min_,'1'))';
% ML=MovieList(allMovieToAnalyse.rawPath(movieIndex(3)),outputPath,'movieListFileName_','metaphaseTest');
% ML=ML.addAnalysisFolder('/project/bioinformatics/Danuser_lab/externBetzig/raw/adavid/lattice/',outputPath)

%%
ML=MovieList.loadMatFile('/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/proudot/u-track-3D/dynROI/testProject/metaphaseTest.mat');


%% Reducing movieSize 
for mIdx=1:ML.getSize()
    MD=ML.getMovie(mIdx);
    MDCrop=crop3D(MD,(MD.getChannel(1).loadStack(1)),'keepFrame',1:5,'name','testDynProjSpeed');
    MD.save();
end


%% Compare running time with project 1D
rehash
for mIdx=1:ML.getSize()
    MD=ML.getMovie(mIdx);
    procs=MD.findProcessTag('crop3D_testDynProjSpeed');
    MDCrop=MovieData.loadMatFile(procs(end).outFilePaths_{1});
    % SpindleRefPack=MD.findPackageName('buildAndProjectSpindleRef');
    pack=MDCrop.searchPackageName('testProjDynROIMetaphasePlate'); 
    testProjDynROIMetaphasePlate(MDCrop,'package',pack(1).setProcess(8,[]));
    MDCrop.save();    
end

%%
MD=ML.getMovie(1);
procs=MD.findProcessTag('crop3D_testDynProjSpeed');
MDCrop=MovieData.loadMatFile(procs(end).outFilePaths_{1});
packs=MDCrop.searchPackageName('testProjDynROIMetaphasePlate');
computeMIPsProcess=packs(1).getProcess(8);

%% Compare GPU running time with project 1D
rehash
for mIdx=1:ML.getSize()
    MD=ML.getMovie(mIdx);
    procs=MD.findProcessTag('crop3D_testDynProjSpeed');
    MDCrop=MovieData.loadMatFile(procs(end).outFilePaths_{1});
    % SpindleRefPack=MD.findPackageName('buildAndProjectSpindleRef');
    pack=MDCrop.searchPackageName('testProjDynROIMetaphasePlate'); 
    %testProjDynROIMetaphasePlate(MDCrop,'package',pack(1).setProcess(4,[]));
    testProjDynROIMetaphasePlate(MDCrop);
    MDCrop.save();    
end