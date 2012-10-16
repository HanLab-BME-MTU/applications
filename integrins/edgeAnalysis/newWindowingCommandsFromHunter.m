%setup movieData object
movieSelectorGUI
movieDataPath = '/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/120828_Cs1C4/analysisCellEdgeSmall3';
MD = MovieData.load(fullfile(movieDataPath,'movieData.mat'));

%determine threshold
thresholdValue = getSegThreshFromFullMovie(MD,1,0.3,1);
close all

%get cell mask
threshParam = struct(...
    'GaussFilterSigma',1,...
    'MethodIndx',1,...
    'ThresholdValue',thresholdValue);
MD = thresholdMovie(MD,threshParam);

%refine cell mask
refinementParam = struct(...
    'ClusureRadius',0,...
    'OpeningRadius',3);
MD = refineMovieMasks(MD,refinementParam);
% MD = refineMovieMasks(MD);

%make movie of mask on top of particle detections
%LOAD MOVIEINFO!
movieMasksParticles(MD,movieInfo,400,1,'movieMasksParticlesOriginal',[],1);

%refine masks using gradient information
threshParamEdge = struct(...
    'filterSigma',1.5,...
    'gradPrctile',[95 90 85 80]);
gapCloseParamEdge = struct(...
    'maxEdgePairDist',13,...
    'maxThetaDiff',pi,...
    'maxC2cAngleThetaDiff',pi/2,...
    'factorContr',[0 0 0 1 0 0],...
    'edgeType',0,...
    'fracImageCell',0.2);
prctileUsed = refineMovieEdgeWithSteerableFilter(MD,threshParamEdge,gapCloseParamEdge,0);
save('paramSteerableFilter','threshParamEdge','gapCloseParamEdge','prctileUsed');

imtool close all

%re-refine for frames where refinement failed
threshParamEdgeRescue = struct(...
    'filterSigma',2.5,...
    'gradPrctile',[95 90 85 80]);
gapCloseParamEdgeRescue = struct(...
    'maxEdgePairDist',13,...
    'maxThetaDiff',pi,...
    'maxC2cAngleThetaDiff',pi/2,...
    'factorContr',[0 0 0 1 0 0],...
    'edgeType',0,...
    'fracImageCell',0.2);
prctileUsedRescue = refineRescueEdgeWithSteerableFilter(MD,threshParamEdgeRescue,gapCloseParamEdgeRescue,1,7);
save('paramSteerableFilterRescue','threshParamEdgeRescue','gapCloseParamEdgeRescue','prctileUsedRescue');

imtool close all

%refine cell mask again with Hunter's code to get back on track
refinementParam2 = struct(...
    'ClosureRadius',0,...
    'OpeningRadius',0);
MD = refineMovieMasks(MD,refinementParam2);

%make movie of mask on top of particle detections
%LOAD MOVIEINFO!
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
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMaskOriginal','AxesHandle',axHandle)
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
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMaskFinal','AxesHandle',axHandle)
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
