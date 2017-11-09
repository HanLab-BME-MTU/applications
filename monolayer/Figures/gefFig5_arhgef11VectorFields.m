function [] = gefFig5_arhgef11VectorFields()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

cntlDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Yunyu_20140911_16hr_5min_0001_AA01_02';
arhgef11Dname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Yunyu_20140911_16hr_5min_0001_AB01_05';

cntlTime = [22];
arhgef11Time = [22];                                              

cntlPatchSize = ceil(15.0/1.25);
arhgef11PatchSize = ceil(15.0/1.25);

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig5/ARHGEF11/';


gefFigUtil_plotVelocityFields(cntlDname,cntlTime,cntlPatchSize,'Cntl',figDname);
gefFigUtil_plotVelocityFields(arhgef11Dname,arhgef11Time,arhgef11PatchSize,'ARHGEF11',figDname);

end


function [] = gefFigUtil_plotVelocityFields(inDname,frames,patchSize,condStr,outDname)

reduceMotionResolution = 0.5;
clim = [0,60];

imgDname = [inDname filesep 'images/']; % xxx.tif 
roiDname = [inDname filesep 'ROI/roi/']; % xxx_roi.mat --> ROI
mfDname = [inDname filesep 'MF/mf/']; % xxx_mf.mat

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
