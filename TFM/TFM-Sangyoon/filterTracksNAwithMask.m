function [tracksNA,trackIdx] = filterTracksNAwithMask(pathForTheMovieDataFile,outputPath,tracksNA,band)
% [tracksNA] = filterTracksNAwithMask(pathForTheMovieDataFile,tracksNA)
% filter out NA tracks, obtain life time of each NA tracks

% input:    pathForTheMovieDataFile:    path to the movieData file (FA
%                   segmentation and NA tracking package should be run beforehand)
% output:   images will be stored in pathForTheMovieDataFile/trackFrames

% Sangyoon Han September 2015

%% Data Set up
% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
MD = MovieData.load(movieDataPath);
% Get whole frame number
nFrames = MD.nFrames_;

nChannels = numel(MD.channels_);
iChan = 0;

outputFilePath = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
try
    load([dataPath filesep 'cropInfo.mat'])
catch
    TFMpackage=MD.getPackage(MD.getPackageIndex('TFMPackage'));
    forceProc =TFMpackage.processes_{4};
    forceField = forceProc.loadChannelOutput;
    cropInfo = [ceil(min(forceField(1).pos(:,1))),ceil(min(forceField(1).pos(:,2))),floor(max(forceField(1).pos(:,1))),floor(max(forceField(1).pos(:,2)))];
end
% Finding which channel has a cell mask information
maskProc = MD.getProcess(MD.getProcessIndex('MaskRefinementProcess'));
for k=1:nChannels
    if maskProc.checkChannelOutput(k)
        iChan = k;
    end
end
% See if there is stage drift correction
iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=MD.processes_{iSDCProc};
    if length(SDCProc.funParams_.ChannelIndex)>1
        iChan = 2;
    elseif length(SDCProc.funParams_.ChannelIndex) == 1
        iChan = SDCProc.funParams_.ChannelIndex;
    else
        error('No channel associated with SDC process!')
    end
    if iChan==2
        iBeadChan=1;
    else
        iBeadChan = SDCProc.funParams_.ChannelIndex(1);
    end
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
end

% Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered
progressText(0,'Filtering adhesion tracks based on cell mask:');
trackIdx = false(numel(tracksNA),1);
for ii=1:nFrames
    % Cell Boundary Mask 
    mask = maskProc.loadChannelOutput(iChan,ii);
    % Apply SDC
    if ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        % Apply subpixel-wise registration to original masks
        I = padarray(mask, [maxY, maxX]);
        maskSDC = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
    else
        maskSDC = mask;
        clear mask
    end
    maskForceGrid = false(size(maskSDC));
    
    maskForceGrid(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band) = ...
        maskSDC(cropInfo(2)+band:cropInfo(4)-band,cropInfo(1)+band:cropInfo(3)-band);
    
    % collect index of tracks which first appear at frame ii
    idxFirstAppear = arrayfun(@(x) x.startingFrame==ii,tracksNA);
    % now see if these tracks ever in the maskOnlyBand
    for k=find(idxFirstAppear)'
        if maskForceGrid(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii)))
            trackIdx(k) = true;
        end
    end
    progressText(ii/nFrames,'Filtering adhesion tracks based on cell mask:');
end    
% get rid of tracks that have out of bands...
tracksNA = tracksNA(trackIdx);