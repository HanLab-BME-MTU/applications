%addpath(genpath('C:' filesep 'Users' filesep 'Philippe' filesep 'code' filesep 'utsw-ssh'))
%load('C:' filesep 'Users' filesep 'Philippe' filesep 'project-local' filesep 'externBetzig' filesep 'analysis' filesep 'adavid' filesep 'smallSample' filesep 'prometaphase' filesep 'earlyCell1_12' filesep 'movieData.mat')
%kinList=[10 32];

load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');
kinList=[2];

%load('C:' filesep 'Users' filesep 'Philippe' filesep 'project-local' filesep 'externBetzig' filesep 'analysis' filesep 'adavid' filesep 'smallSample' filesep 'prometaphase' filesep 'laterCell1_12' filesep 'movieData.mat')
% fixed to pole 1 here

%% loading isotropized trackes in the volume ref.
outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EBTracksRawData=tmp;
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
kinTracksRawData=kinTrackData;

CanonRef=FrameOfRef();
CanonRef.genCanonicalRef(EBTracksRawData.tracksLabRef.numTimePoints());

%% Create Pole and interpolar ref.
% Load addPoleRef to weed out inliers kin
tmp=load(MD.getPackage(1).getProcess(2).outFilePaths_{2});
inliersKin=logical(arrayfun(@(k) k.inliers(1),tmp.kinTracks));
inlierKinTracks=kinTracksRawData.tracksLabRef(inliersKin);

% Load pole detection and build the associated tracks
load(MD.getPackage(1).getProcess(1).outFilePaths_{1});
P1TrackISO=TracksHandle();
P1TrackISO.endFrame=kinTracksRawData.tracksLabRef.numTimePoints;
P1TrackISO.startFrame=1;
P1TrackISO.x=arrayfun(@(d) d.xCoord(1,1),poleMovieInfo)';
P1TrackISO.y=arrayfun(@(d) d.yCoord(1,1),poleMovieInfo)';
P1TrackISO.z=arrayfun(@(d) d.zCoord(1,1)*MD.pixelSizeZ_/MD.pixelSize_,poleMovieInfo)'

P2TrackISO=TracksHandle();
P2TrackISO.endFrame=kinTracksRawData.tracksLabRef.numTimePoints;
P2TrackISO.startFrame=1;
P2TrackISO.x=arrayfun(@(d) d.xCoord(2,1),poleMovieInfo)';
P2TrackISO.y=arrayfun(@(d) d.yCoord(2,1),poleMovieInfo)';
P2TrackISO.z=arrayfun(@(d) d.zCoord(2,1)*MD.pixelSizeZ_/MD.pixelSize_,poleMovieInfo)'

Pole1RefISO=FrameOfRef();
Pole1RefISO.setOriginFromTrack(P1TrackISO);
Pole1RefISO.genCanonicalBase();

InterpolarRefISO=FrameOfRef();
InterpolarRefISO.setOriginFromTrack(P1TrackISO);
InterpolarRefISO.setZFromTrack(P2TrackISO);
InterpolarRefISO.genBaseFromZ();

InterpolarRefISO.applyBaseToTrack(P1TrackISO,'pole1');
InterpolarRefISO.applyBaseToTrack(inlierKinTracks,'pole1');

%% For a few selected kin create frame of reference

for kinIdx=kinList
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'FoF',CanonRef,'crop','full','name',['demoDynManifoldFrameOfRef' filesep 'fullCropLabFrameOfRef-kin' num2str(kinIdx)],'channelRender','grayRed','transType','translation')

%%
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'FoF',CanonRef,'crop','manifold','name',['demoDynManifoldFrameOfRef' filesep 'croppedLabFrameOfRef-kin' num2str(kinIdx)],'channelRender','grayRed','transType','translation')

Pole1RefISO.applyBaseToTrack(P1TrackISO,'P1');
Pole1RefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1');
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1,inlierKinTracks(kinIdx).P1],'FoF',Pole1RefISO,'crop','manifold','name',['demoDynManifoldFrameOfRef' filesep 'croppedPoleFrameOfRef-kin' num2str(kinIdx)],'channelRender','grayRed','transType','translation')
%%

InterpolarRefISO.applyBaseToTrack(P1TrackISO,'P1P2');
InterpolarRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1P2');
project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1P2,inlierKinTracks(kinIdx).P1P2],'FoF',InterpolarRefISO,'crop','manifold','name',['demoDynManifoldFrameOfRef' filesep 'croppedInterpolarFrameOfRef-kin' num2str(kinIdx)],'channelRender','grayRed')

%% estimate
PoleKinRefISO=FrameOfRef();
PoleKinRefISO.setOriginFromTrack(P1TrackISO);
PoleKinRefISO.setZFromTrack(inlierKinTracks(kinIdx));
PoleKinRefISO.genBaseFromZ();
PoleKinRefISO.applyBaseToTrack(P1TrackISO,'P1K');
PoleKinRefISO.applyBaseToTrack(inlierKinTracks(kinIdx),'P1K');

project1D(MD,[P1TrackISO,inlierKinTracks(kinIdx)],'dynPoligonREF',[P1TrackISO.P1K,inlierKinTracks(kinIdx).P1K],'FoF',PoleKinRefISO,'name',['demoDynManifoldFrameOfRef' filesep 'croppedPKinFrameOfRef-kin' num2str(kinIdx)],'crop','manifold','transType','affineOnePass','channelRender','grayRed');
end

%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(2)])
%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(5)])
