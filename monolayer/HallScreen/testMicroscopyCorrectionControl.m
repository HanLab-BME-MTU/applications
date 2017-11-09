function [] = testMicroscopyCorrectionControl()
addpath(genpath('/home2/azaritsky/code/common'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));

dirname = '/project/cellbiology/gdanuser/collab/hall/TestMicroscopyCorrection/';
filename = 'Yunyu_20150501_2hr_5min_0004_AB03_01.zvi';

if ~exist([dirname filename(1:end-4) '/images/'])
    whArrangeData(dirname);
end

params.timePerFrame = 5;
params.nRois = 1;
params.pixelSize = 1.25;

nFrames = 25;
[params,dirs] = whInitParamsDirs(params, dirname, filename(1:end-4), nFrames);

whLocalMotionEstimation(params,dirs);
correctGlobalMotion(params,dirs);

close all;
end

%%
function [] = correctGlobalMotion(params,dirs)
if exist([dirs.mfDataOrig filesep '001_mf.mat'],'file')
    unix(sprintf('cp -R %s %s',[dirs.mfDataOrig '*.mat'],dirs.mfData));
end
unix(sprintf('cp -R %s %s',[dirs.mfData '*.mat'],dirs.mfDataOrig));

correctionsDx = [];
correctionsDy = [];
medianPrecentDx = [];
medianPrecentDy = [];

for t = 1 : params.nTime - params.frameJump
    mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];
    load(mfFname); % dxs, dys
    
    [correctDx, medianPDx] = getCorrection(dxs);
    correctionsDx = [correctionsDx correctDx];
    medianPrecentDx = [medianPrecentDx medianPDx];
    
    [correctDy, medianPDy] = getCorrection(dys);
    correctionsDy = [correctionsDy correctDy];
    medianPrecentDy = [medianPrecentDy medianPDy];
    
    fprintf(sprintf('frame %d: X(%.1f,%.2f), Y(%.1f,%.2f)\n',t,correctDx,medianPDx,correctDx,medianPDx));
end

end

%%
function [correctDx, medianPDx] = getCorrection(dxs)
xsBackground = dxs(:);
xsBackground = xsBackground(~isnan(xsBackground));
nXsBackground = length(xsBackground);
medianXsBackground = median(xsBackground(:));

sumMedian0 = sum(xsBackground == (medianXsBackground-1));
sumMedian1 = sum(xsBackground == medianXsBackground);
sumMedian2 = sum(xsBackground == (medianXsBackground+1));
allSumMedian = sumMedian0 + sumMedian1 + sumMedian2;

correctDx = 0;
if (allSumMedian > 0.6 * nXsBackground)
    correctDx = -(((sumMedian0 * (medianXsBackground-1)) + (sumMedian1 * medianXsBackground) + (sumMedian2 * (medianXsBackground+1)))/allSumMedian);
end

medianPDx = allSumMedian/nXsBackground;
end