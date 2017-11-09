function [] = gefFig1_rhoaVectorFields()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

cntlDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Yunyu_20150401_16hr_5min_0001_AB01_03/';
rhoaDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Yunyu_20150401_16hr_5min_0001_AB03_04/';

cntlTime = [2,14,26];
rhoaTime = [3,13,25]; % 14 and 26 are also options                                               

cntlPatchSize = ceil(15.0/1.25);
rhoaPatchSize = ceil(15.0/1.25);

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig1/velocityFields/';

gefFigUtil_plotVelocityFields(cntlDname,cntlTime,cntlPatchSize,'Cntl',figDname);
gefFigUtil_plotVelocityFields(rhoaDname,rhoaTime,rhoaPatchSize,'Rhoa',figDname);

end


function [] = gefFigUtil_plotVelocityFields(inDname,frames,patchSize,condStr,outDname)

reduceMotionResolution = 0.5;
clim = [0,60];

imgDname = [inDname 'images/']; % xxx.tif 
roiDname = [inDname 'ROI/roi/']; % xxx_roi.mat --> ROI
mfDname = [inDname 'MF/mf/']; % xxx_mf.mat

for itime = 1 : length(frames)
    curTime = frames(itime);
    I = imread([imgDname sprintf('%.3d',curTime) '.tif']); % images data
    load([roiDname sprintf('%.3d',curTime) '_roi.mat']); % ROI
    load([mfDname sprintf('%.3d',curTime) '_mf.mat']); % dxs, dys
    
    outFname = [outDname 'velocityFields_' condStr '_' num2str(itime) '.eps'];
    
    Iinvert = imcomplement(I);
    
    visualizeMotionFieldsColormap(I,ROI,dxs,dys,patchSize,reduceMotionResolution,clim,outFname);
end
end
