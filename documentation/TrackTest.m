% % % connect to Imaris
vImarisApplication = actxserver('Imaris.Application');

% build a surpas scene
vImarisSurpassScene = actxserver('Imaris.DataContainer');
vImarisSurpassScene.AddChild(actxserver('Imaris.LightSource'));
vImarisSurpassScene.AddChild(actxserver('Imaris.Frame'));
vImarisApplication.mSurpassScene = vImarisSurpassScene;

% define some spot coordinates (position um, time index)
vSpotPosXYZ = [ 0.0  0.0  0.0; ... % 0
               10.0 10.0 10.0; ... % 1
               20.0 10.0 15.0; ... % 2
               10.0 20.0 15.0; ... % 3
               10.0 30.0 10.0];    % 4
vSpotPosT = [0;1;2;2;3]; % time index

% create a spots component
vImarisSpots = actxserver('Imaris.Spots');
vImarisSpots.mName = 'my spots';
vImarisSpots.mRadius = 0.3;
vImarisSpots.SetPos(vSpotPosXYZ, vSpotPosT);

% add spots to the surpass scene as last child
vImarisSurpassScene.AddChild(vImarisSpots);

% define edges (build pairs of spots)
vEdges = [0 1;1 2;1 3;3 4;2 4]; % spot index

% create track component
vImarisTrack = actxserver('Imaris.Track');
vImarisTrack.mName = 'my track';
vImarisTrack.SetTrack(vImarisSpots, vEdges);

% add track to the surpass scene as last child
vImarisSurpassScene.AddChild(vImarisTrack);
