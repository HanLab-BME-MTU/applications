function []=reBuildForceMap(MD)
load([MD.movieDataPath_ '/TFMPackage/correctedDisplacementField/displField.mat'])
load([MD.movieDataPath_ '/TFMPackage/forceField/forceField.mat'])
[tMapIn, tmax, tmin, cropInfo,tMapXin,tMapYin] = generateHeatmapShifted(forceField,displField,0);
%% 
movieData=MD;    
nFrames = movieData.nFrames_;

iStep2Proc = movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);
step2Proc = movieData.processes_{iStep2Proc};
pDisp = step2Proc.funParams_;

iProc=4;
forceFieldProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = forceFieldProc.funParams_;
outputFile{1,1} = [p.OutputDirectory filesep 'forceField.mat'];
outputFile{2,1} = [p.OutputDirectory filesep 'tractionMaps.mat'];

%% 
iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
SDCProc=movieData.processes_{iSDCProc};
refFrame = double(imread(SDCProc.outFilePaths_{2,pDisp.ChannelIndex}));
firstMask = refFrame>0; %false(size(refFrame));

disp('Writing traction maps ...')
tMap = cell(1,nFrames);
tMapX = cell(1,nFrames);
tMapY = cell(1,nFrames);
fCfdMap = cell(1,1); %force confidence

% Set up the output directories
outputDir = fullfile(p.OutputDirectory,'tractionMaps');
mkClrDir(outputDir);
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
outFileTMap=@(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
frameSequence=1:nFrames;
tic
for ii=frameSequence
    % starts with original size of beads
    cur_tMap = zeros(size(firstMask));
    cur_tMapX = zeros(size(firstMask));
    cur_tMapY = zeros(size(firstMask));
%     cur_distBeadMap = zeros(size(firstMask));
    cur_tMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{ii};
    cur_tMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapXin{ii};
    cur_tMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapYin{ii};
%     cur_distBeadMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = distBeadMapIn{ii};
    tMap{ii} = cur_tMap;
    tMapX{ii} = cur_tMapX;
    tMapY{ii} = cur_tMapY;
end
toc
disp('Saving ..')
save(outputFile{2},'tMap','tMapX','tMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
disp('Done!')