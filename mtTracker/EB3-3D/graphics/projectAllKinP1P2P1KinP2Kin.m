function projectAllKinP1P2P1KinP2Kin(MD)

relativeOutput='AllKinPrint';
%%
% loading isotropized trackes in the volume ref.
outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EBTracksRawData=tmp;
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
kinTracksRawData=kinTrackData;

% Load addPoleRef to weed out inliers kin
tmp=load(MD.getPackage(1).getProcess(2).outFilePaths_{2});
inliersKin=logical(arrayfun(@(k) k.inliers(1),tmp.kinTracks));
inlierKinTracks=kinTracksRawData.tracksLabRef(inliersKin);

% Load pole detection and build the associated tracks
load(MD.getPackage(1).getProcess(1).outFilePaths_{1});
P1TrackISO=TracksHandle();
P1TrackISO.endFrame=20;
P1TrackISO.startFrame=1;
P1TrackISO.x=arrayfun(@(d) d.xCoord(1,1),poleMovieInfo)';
P1TrackISO.y=arrayfun(@(d) d.yCoord(1,1),poleMovieInfo)';
P1TrackISO.z=arrayfun(@(d) d.zCoord(1,1)*MD.pixelSizeZ_/MD.pixelSize_,poleMovieInfo)'

P2TrackISO=TracksHandle();
P2TrackISO.endFrame=20;
P2TrackISO.startFrame=1;
P2TrackISO.x=arrayfun(@(d) d.xCoord(2,1),poleMovieInfo)';
P2TrackISO.y=arrayfun(@(d) d.yCoord(2,1),poleMovieInfo)';
P2TrackISO.z=arrayfun(@(d) d.zCoord(2,1)*MD.pixelSizeZ_/MD.pixelSize_,poleMovieInfo)'

%% Build the InterpolarRef
InterpolarRefISO=FrameOfRef();
InterpolarRefISO.setOriginFromTrack(P1TrackISO);
InterpolarRefISO.setZFromTrack(P2TrackISO);
InterpolarRefISO.genBaseFromZ();

InterpolarRefISO.applyBaseToTrack(P1TrackISO,'P1P2');
InterpolarRefISO.applyBaseToTrack(P2TrackISO,'P1P2');
InterpolarRefISO.applyBaseToTrack(inlierKinTracks,'P1P2');

for kinIdx=1:length(inlierKinTracks)
InterpolarRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1P2');
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1P2,inlierKinTracks(kinIdx).P1P2],'FoF',InterpolarRefISO,'crop','manifold','name',[relativeOutput filesep 'P1-kin-interpolar-kin-' num2str(kinIdx)],'channelRender','grayRed')
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P2TrackISO.P1P2,inlierKinTracks(kinIdx).P1P2],'FoF',InterpolarRefISO,'crop','manifold','name',[relativeOutput filesep 'P2-kin-interpolar-kin-' num2str(kinIdx)],'channelRender','grayRed')

PoleKinRefISO.setOriginFromTrack(P1TrackISO);
PoleKinRefISO.setZFromTrack(inlierKinTracks(kinIdx));
PoleKinRefISO.genBaseFromZ();
PoleKinRefISO.applyBaseToTrack(P1TrackISO,'P1K');
PoleKinRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1K');
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1K,inlierKinTracks(kinIdx).P1K],'FoF',PoleKinRefISO,'name',[relativeOutput filesep 'P1-kin-static-kin-' num2str(kinIdx)],'channelRender','grayRed')

PoleKinRefISO.setOriginFromTrack(P2TrackISO);
PoleKinRefISO.setZFromTrack(inlierKinTracks(kinIdx));
PoleKinRefISO.genBaseFromZ();
PoleKinRefISO.applyBaseToTrack(P1TrackISO,'P2K');
PoleKinRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P2K');
project1D(MD,[P2TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P2K,inlierKinTracks(kinIdx).P2K],'FoF',PoleKinRefISO,'name',[relativeOutput filesep 'P2-kin-static-kin' num2str(kinIdx)],'channelRender','grayRed')
end
