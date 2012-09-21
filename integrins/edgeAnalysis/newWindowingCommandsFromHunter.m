%setup movieData object
movieSelectorGUI
load movieData.mat

%determine threshold
thresholdValue = getSegThreshFromFullMovie(MD,0.5,0.45,1);
close all

%get cell mask
threshParam.GaussFilterSigma = 0.5;
threshParam.MethodIndx = 1;
threshParam.ThresholdValue = thresholdValue;
MD = thresholdMovie(MD,threshParam);

%refine cell mask
refinementParam.ClusureRadius = 0;
refinementParam.OpeningRadius = 3;
MD = refineMovieMasks(MD,refinementParam);
% MD = refineMovieMasks(MD);

%make movie of mask on top of particle detections
%LOAD MOVIEINFO!
movieMasksParticles(MD,movieInfo,400,1,'movieMasksParticlesOriginal',[],1);

%calculate protrusion vectors
protrusionParam.SegProcessIndex = 2;
MD = getMovieProtrusion(MD,protrusionParam);

%divide mask into windows
windowParam.MethodName = 'ProtrusionBased';
windowParam.ParaSize = 2;
windowParam.PerpSize = 2;
windowParam.SegProcessIndex = 2;
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
