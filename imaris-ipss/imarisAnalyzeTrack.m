function imarisAnalyzeTrack(imarisAppOrID)
% help later


%============
% TEST INPUT
%============

try
    imaApp = imarisCheckHandle(imarisAppOrID);
catch
    error('please specify a valid imaris handle or id')
end

%============



%=======================
% READ DATA FROM IMARIS
%=======================

% find surpass scene
imaSurpassScene = imaApp.mSurpassScene;

% find tracks
numberOfChildren = imaSurpassScene.GetNumberOfChildren;
tracks(1:100) = struct('imaHandle',[],'name',[],...
    'spotXYZ',[], 'spotT',[]);
nTrack = 0;

for kid = 1:numberOfChildren
    % add to list of tracks if we have one
    imaKid = imaSurpassScene.GetChild(kid);
    if imaApp.mFactory.IsTrack(imaKid)
        nTrack = nTrack +1;
        % store handle
        tracks(nTrack).imaHandle = imaKid;
        
        % --- store everything else - we need it, anyway
        tracks(nTrack).name = imaKid.mName;
        
        % all we need is a unique list of coordinates
        edges = imaKid.GetEdges;
        imaSpotHandle = imaKid.GetSpots;
        % convert to double (just in case)
        spotCoord = double(imaSpotHandle.GetPositionsXYZ);
        spotT = double(imaSpotHandle.GetIndicesT);
        
        % edges contains indices to start/end-spot. Remember +1
        connectedSpotIdx = unique(edges(:)+1);
        
        tracks(nTrack).spotXYZ = spotCoord(connectedSpotIdx,:);
        
        % remember +1
        tracks(nTrack).spotT   = spotT(connectedSpotIdx)' + 1;
        
        % --- insert here: control to make sure there is only one spot per 
        %                  timepoint (unique(spotT)) 
        
        clear edges spotCoord spotT connectedSpotIdx
        
    end
end

% remove superfluous entry in track struct
tracks(nTrack+1:end) = [];

% get number of tracks
numTracks = length(tracks);

% if there are no tracks, there's nothing to do here. Die.
if numTracks == 0
    error('no track detected in current Imaris Application')
end


% find and close GUI if exist
h = findall(0,'Tag','iatGUI');
delete(h);

% launch GUI
imarisAnalyzeTrackGUI(tracks);


        