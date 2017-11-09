function [] = runZVIPOC()
pixelSize = 1.249902;%0.624951;
timePerFrame = 5;
initParams.nRois = 1;
mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Hall/POC_Nov21/binning/';
whAnalysis(mainDirname,pixelSize,timePerFrame,initParams);
end