%setup movieData object and fake-analyze to set up the results paths
%correctly
segmentationPackageGUI

%load movieData object
movieDataPath = '/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV/091012_CHO08lofix_mEosAV_2400_25/analysisCellEdgeSmall';
MD = MovieData.load(fullfile(movieDataPath,'movieData.mat'));

%determine threshold
thresholdValue = getSegThreshFromFullMovie(MD,1,0,1);
close all

%get cell mask
threshParam = struct(...
    'GaussFilterSigma',1,...
    'MethodIndx',1,...
    'ThresholdValue',thresholdValue);
MD = thresholdMovie(MD,threshParam);

%refine cell mask
refinementParam = struct(...
    'ClusureRadius',5,...
    'OpeningRadius',3);
MD = refineMovieMasks(MD,refinementParam);

%load detection results
load ../analysisAlphaV/tracks/detection1AllFrames.mat

%make movie of mask on top of particle detections
movieMasksParticles(MD,movieInfo,400,1,'movieMasksParticlesThresh',[],1);

%make images from single molecule signal to enhance edge detection
mkdir imagesSM4Edge
singleMolSignal4Edges('/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/120907/120907_Cs1C2/imagesCellEdge/120907_Cs1C2_CHO_mEos2Av_5minEdgeStack_0001.tif','/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/120907/120907_Cs1C2/imagesAlphaV/120907_Cs1C2_CHO_mEos2Av_00002.tif','/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/120907/120907_Cs1C2/imagesSM4Edge',400,40)

%refine masks using gradient information
threshParamEdge = struct(...
    'filterSigma',1.5,...
    'prctile',[95 90 85 80],...
    'bandHalfWidth',-1);
gapCloseParamEdge = struct(...
    'maxEdgePairDist',9,...
    'factorContr',[1 1 1 1]);
meanBkg = 110;
smDir = '/home/kj35/files/LCCB/receptors/Galbraiths/data/lifeActAndCellEdge/120224/120224_Cs1C1_LifeAct/120224_Cs1C1a_LifeAct/imagesSM4Edge';
prctileUsed = refineMovieEdgeWithSteerableFilter(MD,threshParamEdge,gapCloseParamEdge,1,meanBkg,smDir);
save('paramSteerableFilter','threshParamEdge','gapCloseParamEdge','prctileUsed','meanBkg');

imtool close all

%re-refine for frames where refinement failed
threshParamEdgeRescue = struct(...
    'filterSigma',1.5,...
    'gradPrctile',[90 85 80]);
gapCloseParamEdgeRescue = struct(...
    'maxEdgePairDist',9,...
    'factorContr',[0 1 0 1 1],...
    'edgeType',1,...
    'fracImageCell',0.2);
meanBkg = [];
prctileUsedRescue = refineRescueEdgeWithSteerableFilter(MD,15,threshParamEdgeRescue,gapCloseParamEdgeRescue,0,movieInfo,400,meanBkg);
save('paramSteerableFilterRescue','threshParamEdgeRescue','gapCloseParamEdgeRescue','prctileUsedRescue');

imtool close all

%run mask refinement to get back on track
refinementParam2 = struct(...
    'ClosureRadius',0,...
    'OpeningRadius',0);
MD = refineMovieMasks(MD,refinementParam2);

%make movie of mask on top of particle detections
movieMasksParticles(MD,movieInfo,400,1,'movieMasksParticlesSteerable',[],1);

%fix edges that need fixing manually
hackathonManualSegmentationGUI

%load movieData object again
MD = MovieData.load(fullfile(movieDataPath,'movieData.mat'));

%copy final masks to refined_masks directory
copyHandMasks2Masks(MD)

%run mask refinement to get back on track
refinementParam2 = struct(...
    'ClosureRadius',0,...
    'OpeningRadius',0);
MD = refineMovieMasks(MD,refinementParam2);

%make movie of mask on top of particle detections
movieMasksParticles(MD,movieInfo,400,1,'movieMasksParticlesFinal',[],1);

%calculate protrusion vectors
protrusionParam = struct(...
    'SegProcessIndex',2);
MD = getMovieProtrusion(MD,protrusionParam);

%divide mask into windows
windowParam = struct(...
    'MethodName','ProtrusionBased',...
    'ParaSize',2,...
    'PerpSize',2,...
    'SegProcessIndex',2);
MD = getMovieWindows(MD,windowParam);

%sample protrusion vectors
MD = sampleMovieProtrusion(MD);

%make movie of windows + protrusion vectors
figure('units','normalized','position',[0 0 1 1])
axHandle = gca;
makeMovieMovie(MD,'Overlay','Windows+Protrusion','SegProcessIndex',2,'FileName','movieWindowsProtrusion','AxesHandle',axHandle)




% % %MinMax
% % threshParam.MaxJump = 1.2;
% % threshParam.GaussFilterSigma = 1;
% % threshParam.MethodIndx = 1;
% % MD = thresholdMovie(MD,threshParam);
% % %Otsu
% % threshParam.MaxJump = 1.2;
% % threshParam.GaussFilterSigma = 1;
% % threshParam.MethodIndx = 2;
% % MD = thresholdMovie(MD,threshParam);
% % %Rosin
% % threshParam.MaxJump = 1.2;
% % threshParam.GaussFilterSigma = 1;
% % threshParam.MethodIndx = 3;
% % MD = thresholdMovie(MD,threshParam);
% % %Gradient
% % threshParam.MaxJump = [];
% % threshParam.GaussFilterSigma = [];
% % threshParam.MethodIndx = 4;
% % MD = thresholdMovie(MD,threshParam);
% %Fixed threshold
% % threshParam.MaxJump = [];
% 
%make movie of mask
figure('units','normalized','position',[0 0 1 1])
axHandle = gca;
makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMaskOriginal','AxesHandle',axHandle)
% 
% %refine masks using gradient information
% closureRadius = 5;
% edgeThresh = refineMovieEdgeWithSteerableFilter(MD,1,closureRadius);
% save('paramSteerableFilter','closureRadius','edgeThresh');
% 
% imtool close all
% 
% %refine cell mask again with Hunter's code to get back on track
% refinementParam2.ClosureRadius = 0;
% refinementParam2.OpeningRadius = 0;
% MD = refineMovieMasks(MD,refinementParam2);
% 
% %make movie of final mask
figure('units','normalized','position',[0 0 1 1])
axHandle = gca;
makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMaskFinal','AxesHandle',axHandle)
% 
% %make movie of mask on top of particle detections
% %LOAD MOVIEINFO!
% movieMasksParticles(movieInfo,400,[],[],1,'movieMasksParticlesFinal',MD.movieDataPath_,[],1);
% 
% %make movie of protrusion vectors
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Protrusion','SegProcessIndex',2,'FileName','movieProtrusion','AxesHandle',axHandle)
% 
% %make movie of windows
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Windows','SegProcessIndex',2,'FileName','movieWindows','AxesHandle',axHandle)
% 
% % refinementParam.MaxEdgeGap = 20;
% % refinementParam.PreEdgeGrow = 0;
% % refinementParam.MaxEdgeAdjust = 20;
% % refinementParam.EdgeRefinement = 1;
