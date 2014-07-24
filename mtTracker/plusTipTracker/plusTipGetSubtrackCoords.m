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
%
% Kathryn Applegate, 2010
% Sebastien Besson, last modified Feb 2012

if nargin<3 || isempty(useMerged) || useMerged~=1
    trackData=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix; % all
else
    trackData=plusTipMergeSubtracks(projData); % first output is merged
end

% get data from all subtracks
if nargin>=2 && ~isempty(idx), trackData=trackData(idx,:); end

% total number of subtracks
nSubtracks=size(trackData,1);
% list of corresponding track numbers for each subtrack
trackNum=trackData(:,1);
% list of start frame for each subtrack
sF=trackData(:,2);
% list of end frame for each subtrack
eF=trackData(:,3);
% Create linear frame index 
frameIndex=[0; cumsum(eF-sF+1)];

% initialize where coordinates will be stored
xMat=NaN(nSubtracks,projData.nFrames);
yMat=NaN(nSubtracks,projData.nFrames);

% for each subtrack, get list of frames over which subtrack exists
framesPerSubtrack=zeros(frameIndex(end),1);
for i=1:nSubtracks, framesPerSubtrack(frameIndex(i)+1:frameIndex(i+1))=sF(i):eF(i); end

% for each frame of each subtrack, write corresponding TRACK number
trackNumPerSub=zeros(frameIndex(end),1);
for i=1:nSubtracks, trackNumPerSub(frameIndex(i)+1:frameIndex(i+1))=trackNum(i); end
idx=sub2ind([size(projData.xCoord,1), projData.nFrames],trackNumPerSub,framesPerSubtrack);

% for each frame of each subtrack, write corresponding SUBTRACK number
subNumPerSub=zeros(frameIndex(end),1);
for i=1:nSubtracks, subNumPerSub(frameIndex(i)+1:frameIndex(i+1))=i; end
idx2=sub2ind(size(xMat),subNumPerSub,framesPerSubtrack);

% Fill the matrix
xMat(idx2)=projData.xCoord(idx);
yMat(idx2)=projData.yCoord(idx);