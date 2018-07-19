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
ip.addParameter('InputMicrons',false);%If true, tracks are input with coordinates in microns. If not will be adjusted for display
ip.addParameter('Name','Tracks');%What to name the track objects
ip.addParameter('ShowEdges',true);%Show tracked edges (why would I add this feature? because I'm in a rush and this function sucks anyways.... This is a workaround for missing detections
ip.parse(varargin{:});
p = ip.Results;


% -- setup the spots --- %

spotHandle = iceConn.mImarisApplication.GetFactory.CreateSpots;

nPtsPer = arrayfun(@(t)numel(t.x),tracks);

%All spot coordinates, times
voxSize = iceConn.getVoxelSizes;
X = vertcat([tracks.x],[tracks.y],[tracks.z])';
if ~p.InputMicrons
    X = X * voxSize(1);%The voxel anisotropy has been corrected in these detections, so we use only XY pixel size.
end

T = [tracks.t]';
R = ones(sum(nPtsPer),1)*p.SpotRadius;
nanX = any(isnan(X),2);
spotHandle.Set(X(~nanX,:),T(~nanX),R(~nanX));

% --- setup the track edges ---- %
if p.ShowEdges
    nTracks = numel(tracks);
    trackEdges = cell(nTracks,1);
    startInd = 0;
    for j = 1:nTracks

        trackEdges{j} = [(startInd:startInd+nPtsPer(j)-2)' (startInd+1:startInd+nPtsPer(j)-1)'];
        startInd = startInd + nPtsPer(j);

    end

    spotHandle.SetTrackEdges(vertcat(trackEdges{:}));
end
% -- add to imaris scene --- %

% Set the name
spotHandle.SetName(p.Name);
%add to scene
iceConn.mImarisApplication.GetSurpassScene.AddChild(spotHandle, -1);