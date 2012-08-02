% %setup movieData object
% movieSelectorGUI
% load movieData.mat
% 
% %get cell mask
% %Rosin
% threshParam.MaxJump = 1.2;
% threshParam.GaussFilterSigma = 0.5;
% threshParam.MethodIndx = 3;
% MD = thresholdMovie(MD,threshParam);
% %MinMax
% threshParam.MaxJump = 1.2;
% threshParam.GaussFilterSigma = 0.5;
% threshParam.MethodIndx = 1;
% MD = thresholdMovie(MD,threshParam);
% 
% %refine cell mask
% % refinementParam.MaxEdgeAdjust = 50;
% % refinementParam.MaxEdgeGap = 20;
% % refinementParam.PreEdgeGrow = 0;
% % refinementParam.EdgeRefinement = 1;
% % refinementParam.EdgeRefinement = 0;
% % refinementParam.ClosureRadius = 0;
% % MD = refineMovieMasks(MD,refinementParam);
% MD = refineMovieMasks(MD);
% 
% %make movie of mask
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMaskOriginal','AxesHandle',axHandle)
% 
% %refine masks using gradient information
% refineMovieEdgeWithSteerableFilter(MD,1)
% 
% imtool close all
% 
% %refine cell mask again with Hunter's code to get back on track
% % refinementParam.ClosureRadius = 5;
% % MD = refineMovieMasks(MD,refinementParam);
% MD = refineMovieMasks(MD);
% 
% %make movie of final mask
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Mask','SegProcessIndex',2,'FileName','movieMaskFinal','AxesHandle',axHandle)
% 
% %calculate protrusion vectors
% protrusionParam.SegProcessIndex = 2;
% MD = getMovieProtrusion(MD,protrusionParam);
% 
% %make movie of protrusion vectors
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Protrusion','SegProcessIndex',2,'FileName','movieProtrusion','AxesHandle',axHandle)
% 
% %divide mask into windows
% tic
windowParam.MethodName = 'ProtrusionBased';
windowParam.ParaSize = 2;
windowParam.PerpSize = 2;
windowParam.SegProcessIndex = 2;
MD = getMovieWindows(MD,windowParam);
% toc
% 
% %make movie of windows
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Windows','SegProcessIndex',2,'FileName','movieWindows','AxesHandle',axHandle)
% 
% %make movie of windows + protrusion vectors
% tic
% figure('units','normalized','position',[0 0 1 1])
% axHandle = gca;
% makeMovieMovie(MD,'Overlay','Windows+Protrusion','SegProcessIndex',2,'FileName','movieWindowsProtrusion','AxesHandle',axHandle)
% toc

%sample protrusion vectors
MD = sampleMovieProtrusion(MD);
