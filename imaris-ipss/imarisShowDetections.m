function spotHandle = imarisShowDetections(movieInfo,iceConn,varargin)
%IMARISSHOWDETECTIONS shows detected particles in imaris
%
% spotHandles = imarisShowDetections(movieInfo,iceConn);
%
% Detections should be input in the movieinfo structure format (as used in
% e.g. uTrack)
%
%Hunter Elliott
%11/2016

ip = inputParser;
ip.addParameter('SpotRadius',.25);%Radius of spots indicating detections
ip.parse(varargin{:});
p = ip.Results;


% -- setup the spots --- %

spotHandle = iceConn.mImarisApplication.GetFactory.CreateSpots;

nPtsPer = arrayfun(@(t)size(t.xCoord,1),movieInfo);%Num detections per frame

%All spot coordinates, times
voxSize = iceConn.getVoxelSizes;
X = [vertcat(movieInfo.xCoord) vertcat(movieInfo.yCoord) vertcat(movieInfo.zCoord)]; 
X = X(:,1:2:end);%We don't use positional uncertainties
X = bsxfun(@times,X,voxSize);%Adjust for pixel size

T = arrayfun(@(x)(x*ones(nPtsPer(x),1)),1:numel(movieInfo),'Unif',0);
T = vertcat(T{:})-1;%Adjust for zero indexing
R = ones(sum(nPtsPer),1)*p.SpotRadius;
spotHandle.Set(X,T,R);


% -- add to imaris scene --- %

% Set the name
spotHandle.SetName('Detections');
%Add to scene so it's visible
iceConn.mImarisApplication.GetSurpassScene.AddChild(spotHandle, -1);