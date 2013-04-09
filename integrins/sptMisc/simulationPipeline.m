
% topDir = 'C:\kjData\Galbraiths\data\simulations\mimicFarn\exampleProt2Retr0p5_05';
topDir = 'C:\kjData\Galbraiths\data\simulations\mimicAlphaV\exampleProt2Retr2_05';
topDirF = 'C:\kjData\Galbraiths\data\simulations\mimicFarn\exampleProt2Retr2_05';
cd(topDir)

%make necessary directories
% mkdir analysisFarnSim
mkdir analysisAlphaVSim
mkdir analysisCellEdgeSmall
% mkdir imagesCellEdge
% mkdir imagesCellEdgeFull

% %generate cell images
% cd imagesCellEdgeFull
% dir2save = [topDir '\imagesCellEdgeFull'];
% fileName = 'imagesCellEdgeSim';
% pos0 = 14;
% totalTime = 200;
% sampInt = 0.025;
% pixSize = 0.111;
% imSize = [256 256];
% modelParam = struct('protSpeed',2,'retrSpeed',0.5,'protTime',[100 1],'retrTime',[50 1]);
% simCellEdgeSimple(modelParam,pos0,totalTime,sampInt,pixSize,imSize,dir2save,fileName);
% save('simParamEdge','pos0','totalTime','sampInt','pixSize','imSize','modelParam');
% 
% %sub-sample cell images like in the experimental data
% sourceDir = [topDir '\imagesCellEdgeFull'];
% targetDir = [topDir '\imagesCellEdge'];
% copySubsetOfFiles({sourceDir,targetDir},1:400:8001,'imagesCellEdgeSim_','tif',5,1);

%generate tracks
% cd ../analysisFarnSim
cd analysisAlphaVSim
mkdir tracks
mkdir diffusion
mkdir furtherAnalysis
cd tracks
imSize = [256 256];
numP = 220;
numF = 8000;
intVec = [1 0.1];
load ../../../../lftDistr.mat
% motionParam = struct('diffCoef2D',[1 1.2],'confRad2D',[10 10],'speed1D',[0 0],'probVelSwitch',0,'fracType',[0 1 0]);
motionParam = struct('diffCoef2D',[0.3 0.5],'confRad2D',[10 10],'speed1D',[0 0],'probVelSwitch',0,'fracType',[0 1 0]);
[~,tracksSim] = simulateRandomPlusLinearMotion(imSize,numP,lftDistr,numF,intVec,motionParam);
save('simParamTracks','imSize','numP','numF','intVec','motionParam');
save('tracksSim0','tracksSim');

%retain only tracks in cell at each frame
% firstMaskFile = fullfile(topDir,'imagesCellEdgeFull','imagesCellEdgeSim_00001.tif');
firstMaskFile = fullfile(topDirF,'imagesCellEdgeFull','imagesCellEdgeSim_00001.tif');
indxTracksInCellMask = findTracksInCellMask(tracksSim,firstMaskFile,1:8001);
tracksFinal = tracksSim(indxTracksInCellMask);
save('tracks1AllFrames','tracksFinal','indxTracksInCellMask')

%make fake diffusion analysis
cd ../diffusion
load ../../../../diffAnalysisRes1.mat
diffAnalysisRes = repmat(diffAnalysisRes1,length(tracksFinal),1);
save('diffAnalysis1AllFrames','diffAnalysisRes');

%analyze cell edge activity as for experimental data

%analyze single particle behavior as with experimental data, starting with
%retaining tracks in cell mask, etc.