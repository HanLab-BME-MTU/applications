addpath(genpath('C:\Users\Philippe\code\utsw-ssh'))
load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat')
%kinIdx=10;
kinIdx=2;

%load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\laterCell1_12\movieData.mat')
% fixed to pole 1 here

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

P1TrackANI=P1TrackISO.copy();
P1TrackANI.z=P1TrackANI.z*MD.pixelSize_/MD.pixelSizeZ_;

P2TrackANI=P2TrackISO.copy();
P2TrackANI.z=P2TrackANI.z*MD.pixelSize_/MD.pixelSizeZ_;

Pole1Ref=FrameOfRef();
Pole1Ref.setOriginFromTrack(P1TrackANI);
Pole1Ref.genCanonicalBase();

Pole1RefISO=FrameOfRef();
Pole1RefISO.setOriginFromTrack(P1TrackISO);
Pole1RefISO.genCanonicalBase();

InterpolarRef=FrameOfRef();
InterpolarRef.setOriginFromTrack(P1TrackANI);
InterpolarRef.setZFromTrack(P2TrackANI);
InterpolarRef.genBaseFromZ();


InterpolarRefISO=FrameOfRef();
InterpolarRefISO.setOriginFromTrack(P1TrackISO);
InterpolarRefISO.setZFromTrack(P2TrackISO);
InterpolarRefISO.genBaseFromZ();

InterpolarRefISO.applyBaseToTrack(P1TrackISO,'pole1');
InterpolarRefISO.applyBaseToTrack(inlierKinTracks,'pole1');



%%
%project1D(MD,[P1TrackISO,inlierKinTracks(10)],'crop','full','name','fullKinI10')
% project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'crop','manifold','name','manifKinI10')

%% testing single pass
kinIdx=32;
kinIdx=10;

%%
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'FoF',CanonRef,'crop','full','name',['onePassTesting\fullCropLabFOFAffineOnePass-kin' num2str(kinIdx)],'channelRender','grayRed','transType','affineOnePass')
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'FoF',CanonRef,'crop','manifold','name',['onePassTesting\croppedLabFoFAffineOnePass-kin' num2str(kinIdx)],'channelRender','grayRed','transType','affineOnePass')
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'FoF',CanonRef,'crop','full','name',['onePassTesting\fullCropLabFOFTrans-kin' num2str(kinIdx)],'channelRender','grayRed','transType','translation')
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'FoF',CanonRef,'crop','manifold','name',['onePassTesting\croppedLabFOFTrans-kin' num2str(kinIdx)],'channelRender','grayRed','transType','translation')
%%


Pole1RefISO.applyBaseToTrack(P1TrackISO,'P1');
Pole1RefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1');
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1,inlierKinTracks(kinIdx).P1],'FoF',Pole1RefISO,'crop','manifold','name',['onePassTesting\croppedPoleFrameOfRef-kin' num2str(kinIdx)],'channelRender','grayRed','transType','translation')
%%
InterpolarRefISO.applyBaseToTrack(P1TrackISO,'P1P2');
InterpolarRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1P2');
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1P2,inlierKinTracks(kinIdx).P1P2],'FoF',InterpolarRefISO,'crop','manifold','name',['onePassTesting\croppedInterpolarFrameOfRef-kin' num2str(kinIdx)],'channelRender','grayRed')
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1P2,inlierKinTracks(kinIdx).P1P2],'FoF',InterpolarRefISO,'crop','manifold','name',['onePassTesting\fullInterpolarFrameOfRefAffineOnePass-kin' num2str(kinIdx)],'channelRender','grayRed','transType','affineOnePass')


%% Manifold referential 
PoleKinRefISO=FrameOfRef();
PoleKinRefISO.setOriginFromTrack(P1TrackISO);
PoleKinRefISO.setZFromTrack(inlierKinTracks(kinIdx));
PoleKinRefISO.genBaseFromZ();
PoleKinRefISO.applyBaseToTrack(P1TrackISO,'P1K');
PoleKinRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1K');
%%
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1K,inlierKinTracks(kinIdx).P1K],'FoF',PoleKinRefISO,'name',['onePassTesting\croppedPKinFrameOfRef-kin' num2str(kinIdx)],'crop','manifold','tracks',PoleKinRefISO.applyBase(EBTracksRawData.tracksLabRef,'P1KIn'),'transType','affineOnePass','channelRender','grayRed');


%%
%kinIdx=36
kinIdx=32;
for kinIdx=2:2:40
PoleKinRefISO=FrameOfRef();
PoleKinRefISO.setOriginFromTrack(P1TrackISO);
PoleKinRefISO.setZFromTrack(inlierKinTracks(kinIdx));
PoleKinRefISO.genBaseFromZ();

PoleKinRefISO.applyBaseToTrack(P1TrackISO,'P1K');
PoleKinRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1K');

project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1K,inlierKinTracks(kinIdx).P1K],'FoF',PoleKinRefISO,'name',['P1KinRefAffine2p-kin-' num2str(kinIdx)],'crop','manifold','transType','affine');
end

%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(2)])
%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(5)])
