function [] = runHallPOC()
% pixelSize = 0.624951;
timePerFrame = 5;
initParams.nRois = 1;
mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Hall/POC_Nov21/';
load([mainDirname 'GTPasesScreenMetaData.mat']); % exps
whAnalysis(mainDirname,exps,timePerFrame,initParams);
end