function [] = testScript(filename,params)

addpath(genpath('/project/bioinformatics/Danuser_lab/GEFscreen/analysis/testSW/')); % this is the folder where the code is at.
addpath(genpath('/home2/azaritsky/code/extern'));
addpath(genpath('/home2/azaritsky/code/common'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));

if nargin < 1
    filename = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/testSW/Angeles_20150402_14hrs_5min_AA01_7.tif';
end

if nargin < 2
    params.timePerFrame = 5;
    params.nRois = 1;
    params.always = false;
    params.pixelSize = 1;
    params.patchSize = 15.0; % um
    params.nTime = floor(200 / params.timePerFrame); % frames to process (200 minutes)
    params.maxSpeed = 90; % um / hr (max cell speed)
    
    % for kymographs display
    params.kymoResolution.maxDistMu = 180; % um
    params.kymoResolution.min = params.patchSize;
    params.kymoResolution.stripSize = params.patchSize;
    params.kymoResolution.max = ceil(params.kymoResolution.maxDistMu/params.pixelSize);     
end

processTimeLapse(filename, params);
close all;
end