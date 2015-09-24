function [tracksNA] = filterTracksNAwithMask(pathForTheMovieDataFile,tracksNA)
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


bandArea = zeros(nFrames,1);
nChannels = numel(MD.channels_);
iChan = 0;
% Finding which channel has a cell mask information
maskProc = MD.getProcess(MD.getProcessIndex('MaskRefinementProcess'));
for k=1:nChannels
    if maskProc.checkChannelOutput(k)
        iChan = k;
    end
end
% Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered
disp(['Filtering adhesion tracks based on cell mask. Adhesions appearing at the edge are only considered'])
trackIdx = false(numel(tracksNA),1);
bandwidthNA = 7; %um 
bandwidthNA_pix = round(bandwidthNA*1000/MD.pixelSize_);
for ii=1:nFrames
    % Cell Boundary Mask 
    mask = maskProc.loadChannelOutput(iChan,ii);
    % mask for band from edge
    iMask = imcomplement(mask);
    distFromEdge = bwdist(iMask);
    bandMask = distFromEdge <= bandwidthNA_pix;

    maskOnlyBand = bandMask & mask;
    bandArea(ii) = sum(maskOnlyBand(:)); % in pixel

    % collect index of tracks which first appear at frame ii
    idxFirstAppear = arrayfun(@(x) x.startingFrame==ii,tracksNA);
    % now see if these tracks ever in the maskOnlyBand
    for k=find(idxFirstAppear)'
        if maskOnlyBand(round(tracksNA(k).yCoord(ii)),round(tracksNA(k).xCoord(ii)))
            trackIdx(k) = true;
        end
    end
end    
% get rid of tracks that have out of bands...
tracksNA = tracksNA(trackIdx);