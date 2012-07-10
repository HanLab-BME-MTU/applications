clear all
cd tracks
load tracks1AllFrames.mat
cd ../diffusion
load diffAnalysis1AllFrames.mat
[fracTrajClass,f2fDispRMS] = plotPropertySpatialMap2D(tracksFinal,...
    '120608_Cs3C2_Av717',5,[1 2 3 5],1,[],[],diffAnalysisRes);
save('summary1','f2fDispRMS','fracTrajClass');
