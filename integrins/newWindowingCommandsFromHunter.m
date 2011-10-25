%setup movieData object
movieSelectorGUI
load movieData.mat

%get cell mask
threshParam.ThresholdValue = 570;
% threshParam.GaussFilterSigma = 3;
% threshParam.MaxJump = 0.2;
MD = thresholdMovie(MD,threshParam);
% MD = thresholdMovie(MD);

%refine cell mask
refinementParam.MaxEdgeAdjust = 50;
% refinementParam.MaxEdgeGap = 20;
refinementParam.PreEdgeGrow = 0;
refinementParam.EdgeRefinement = 1;
% refinementParam.EdgeRefinement = 0;
% refinementParam.ClosureRadius = 3;
MD = refineMovieMasks(MD,refinementParam);
% MD = refineMovieMasks(MD);

%make movie of mask
figure('units','normalized','position',[0 0 1 1])
axHandle = gca;
makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMask','AxesHandle',axHandle)

%calculate protrusion vectors
protrusionParam.SegProcessIndex = 2;
MD = getMovieProtrusion(MD,protrusionParam);

%make movie of protrusion vectors
figure('units','normalized','position',[0 0 1 1])
axHandle = gca;
makeMovieMovie(MD,'Overlay','Protrusion','SegProcessIndex',2,'FileName','movieProtrusion','AxesHandle',axHandle)

%divide mask into windows
windowParam.MethodName = 'ProtrusionBased';
windowParam.ParaSize = 5;
windowParam.PerpSize = 5;
windowParam.SegProcessIndex = 2;
MD = getMovieWindows(MD,windowParam);

%make movie of windows
figure('units','normalized','position',[0 0 1 1])
axHandle = gca;
makeMovieMovie(MD,'Overlay','Windows','SegProcessIndex',2,'FileName','movieWindows','AxesHandle',axHandle)

%sample protrusion vectors
MD = sampleMovieProtrusion(MD);
