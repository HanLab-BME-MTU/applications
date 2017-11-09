function [] = runErik1m_stable()
pixelSize = 0.325;
timePerFrame = 1;
mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/131122/2/';
pcAnalysis_stable(mainDirname,pixelSize,timePerFrame);
end