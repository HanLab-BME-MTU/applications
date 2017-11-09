function [] = gefFig5_arhgef18VectorFields()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

cntlDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Yunyu_20130918_20hr_5min_2013_0001_pS1_02';
arhgef18Dname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Yunyu_20130918_20hr_5min_2013_0001_172_02';

% cntlTime = [1,13,25];
cntlTime = [2,13,25];
arhgef18Time = [2,14,27]; % 14 and 26 are also options                                               

cntlPatchSize = ceil(15.0/0.625);
arhgef18PatchSize = ceil(15.0/0.625);

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig5/ARHGEF18/';


warning('something is wrong with the control vector fields');

gefFigUtil_plotVelocityFields(cntlDname,cntlTime,cntlPatchSize,'Cntl',figDname);
gefFigUtil_plotVelocityFields(arhgef18Dname,arhgef18Time,arhgef18PatchSize,'ARHGEF18',figDname);

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
