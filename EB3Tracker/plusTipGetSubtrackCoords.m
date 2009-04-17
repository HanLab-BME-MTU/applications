function [xMat,yMat]=plusTipGetSubtrackCoords(projData,idx)
% return nIdx x nFrames matrices containing xy-coordinates for subtracks
% coords in frames where subtracks do not exist are backfilled with NaNs

% get data from all subtracks
if nargin<2 || isempty(idx)
    trackData=projData.nTrack_start_end_velMicPerMin_class_lifetime; % all
else
    trackData=projData.nTrack_start_end_velMicPerMin_class_lifetime(idx,:);
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


