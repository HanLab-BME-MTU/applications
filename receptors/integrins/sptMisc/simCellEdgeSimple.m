function simCellEdgeSimple(modelParam,pos0,totalTime,sampInt,pixSize,imSize,...
    dir2save,fileName)
%SIMCELLEDGESIMPLE simulates the protrusion and retraction of a simple linear edge
%
%SYNOPSIS simCellEdgeSimple(protSpeed,retrSpeed,protTime,retrTime)
%
%INPUT  modelParam  : Structure containing model parameters:
%           .protSpeed: Protrusion speed [um/min].
%           .retrSpeed: Retraction speed [um/min].
%           .protTime : 2-element row vector containing mean and
%                       standard deviation of protrusion time [s].
%           .retrTime : 2-element row vector containing mean and
%                       standard deviation of shrinkage time [s].
%       pos0        : Initial position of cell edge in image [um].
%       totalTime   : Total time of simulation [s].
%       sampInt     : Sampling interval, i.e. time between frames [s].
%       pixSize     : Pixel size in image [um]
%       imSize      : Image size in X and Y [pixels].
%       dir2save    : Directory where images are to be saved.
%       fileName    : Name of image files. If name is X, image files will
%                     be X_00001, X_00002, etc.
%       
%OUTPUT A series of images indicating cell edge location.
%
%REMARKS Simulation assumes a gamma distribution for protrusion and retraction time.
%        In fact the code calls mtGammaTimeDistr to simulate edge movement.
%
%Khuloud Jaqaman, February 2013

%% Edge position simulation

%map nomenclature
%make protrusion = shrinkage and retraction = growth so that the cell edge
%always starts with a retraction - this inversion will be corrected later
modelParamMT.growthSpeed = modelParam.retrSpeed;
modelParamMT.shrinkageSpeed = modelParam.protSpeed;
modelParamMT.growthTime = modelParam.retrTime;
modelParamMT.shrinkageTime = modelParam.protTime;

%generate edge trajectory
simTraj = mtGammaTimeDistr(modelParamMT,imSize(1)*pixSize-pos0,totalTime);
simTraj(:,2) = imSize(1)*pixSize - simTraj(:,2);

%sample trajectory at proper frames
sampTraj = sampleTraj(simTraj,sampInt);

%retain only frames of interest
numFrames = ceil(totalTime/sampInt)+1;
sampTraj = sampTraj(1:numFrames,:);

%convert edge position to pixels
edgePos = ceil(sampTraj(:,2)/pixSize);

%display trajectory
plot(sampTraj(:,1),edgePos);

%% Image generation

%go over each time point and make images
for iFrame = 1 : numFrames
    
    %initialize image
    image = zeros(imSize);
    
    %add cell to image
    image(1:edgePos(iFrame),:) = 1;
    
    %save image
    imwrite(image,fullfile(dir2save,[fileName '_' num2str(iFrame,'%05i') '.tif']),'tif');
    
end