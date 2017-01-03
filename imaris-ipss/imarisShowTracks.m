function spotHandle = imarisShowTracks(tracks,iceConn,varargin)
%IMARISSHOWTRACKS shows tracked particles and tracks in imaris
%
% spotHandles = imarisShowTracks(tracks,iceConn);
%
% Tracks should be input as a TracksHandle object
%
%Hunter Elliott
%11/2016

ip = inputParser;
ip.addParameter('SpotRadius',.25);%Radius of spots indicating detections
ip.parse(varargin{:});
p = ip.Results;


% -- setup the spots --- %

spotHandle = iceConn.mImarisApplication.GetFactory.CreateSpots;

nPtsPer = arrayfun(@(t)numel(t.x),tracks);

%All spot coordinates, times
voxSize = iceConn.getVoxelSizes;
X = vertcat([tracks.x],[tracks.y],[tracks.z])' * voxSize(1);%The voxel anisotropy has been corrected in these detections, so we use only XY pixel size.

T = [tracks.t]';
R = ones(sum(nPtsPer),1)*p.SpotRadius;
spotHandle.Set(X,T,R);

% --- setup the track edges ---- %
nTracks = numel(tracks);
trackEdges = cell(nTracks,1);
startInd = 0;
for j = 1:nTracks
    
    trackEdges{j} = [(startInd:startInd+nPtsPer(j)-2)' (startInd+1:startInd+nPtsPer(j)-1)'];
    startInd = startInd + nPtsPer(j);
    
end

spotHandle.SetTrackEdges(vertcat(trackEdges{:}));

% -- add to imaris scene --- %

% Set the name
spotHandle.SetName('Tracks');
%add to scene
iceConn.mImarisApplication.GetSurpassScene.AddChild(spotHandle, -1);