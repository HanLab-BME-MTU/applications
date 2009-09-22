function [xMat,yMat]=plusTipGetSubtrackCoords(projData,idx,useMerged)
% plusTipGetSubtrackCoords returns subtrack coordinates
%
% INPUT : projData : project data from meta folder
%         idx      : vector containing subtrack indices for which
%                    coordinates are needed
%         useMerged: 1 if the indices should correspond to the data matrix
%                    where unconventional track types are consolidated into
%                    growth events (e.g. trackTypes=5 are merged with
%                    flanking growth phases)
% OUTPUT: xMat,yMat: nIndex x nFrames matrices containing track coordinates
%                    in frames where subtracks do not exist, matrices are
%                    backfilled with NaNs
%
% this function is called by:
% plusTipParamPlot
% plusTipEventSpeedOverlay
% plusTipPlotRandTraj
% plusTipPlotTracks


if nargin<3 || isempty(useMerged) || useMerged~=1
    trackData=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix; % all
else
    trackData=plusTipMergeSubtracks(projData); % first output is merged
end

% get data from all subtracks
if nargin<2 || isempty(idx)
    % we already have all the data
else
    % just take the subset from idx
    trackData=trackData(idx,:);
end

% list of corresponding track numbers for each subtrack
trackNum=trackData(:,1);
% list of start frame for each subtrack
sF=trackData(:,2);
% list of end frame for each subtrack
eF=trackData(:,3);
% total number of subtracks
nSubtracks=size(trackData,1);

% initialize where coordinates will be stored
xMat=nan(nSubtracks,projData.numFrames);
yMat=nan(nSubtracks,projData.numFrames);

% for each subtrack, get list of frames over which subtrack exists
framesPerSubtrack=arrayfun(@(x,y) [x:y]', sF,eF,'UniformOutput',0);
% for each subtrack, get length in frames
len=cell2mat(arrayfun(@(x) length(x{:}), framesPerSubtrack,'UniformOutput',0));

% for each frame of each subtrack, write corresponding TRACK number
trackNumPerSub=arrayfun(@(x,y) x*ones(1,y)', trackNum,len,'UniformOutput',0);
% convert i,j to index
idx=cellfun(@(x,y) sub2ind([projData.numTracks, projData.numFrames],x,y), trackNumPerSub,framesPerSubtrack,'UniformOutput',0);
% coordinates for features in all subtracks at all frames
coordsX=cellfun(@(i) projData.xCoord(i), idx,'UniformOutput',0);
coordsY=cellfun(@(i) projData.yCoord(i), idx,'UniformOutput',0);


% for each frame of each subtrack, write corresponding SUBTRACK number
subNumPerSub=arrayfun(@(x,y) x*ones(1,y)', [1:nSubtracks]',len,'UniformOutput',0);
% convert i,j to index
idx=cellfun(@(x,y) sub2ind(size(xMat),x,y),subNumPerSub,framesPerSubtrack,'UniformOutput',0);
% fill in matrices with coordinates
xMat(cell2mat(idx))=cell2mat(coordsX);
yMat(cell2mat(idx))=cell2mat(coordsY);

