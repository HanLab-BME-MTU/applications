addpath(genpath('C:\Users\Philippe\code\utsw-ssh'))
load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat')

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

%%
InterpolarRefISO=FrameOfRef();
InterpolarRefISO.setOriginFromTrack(P1TrackISO);
InterpolarRefISO.setZFromTrack(P2TrackISO);
InterpolarRefISO.genBaseFromZ();

InterpolarRefISO.applyBaseToTrack(P1TrackISO,'pole1');
InterpolarRefISO.applyBaseToTrack(inlierKinTracks,'pole1');

%%
%project1D(MD,[P1TrackISO,inlierKinTracks(10)],'crop','full','name','fullKinI10')
project1D(MD,[P1TrackISO,inlierKinTracks(10)],'crop','manifold','name','manifKinI10')
project1D(MD,[P1TrackISO,inlierKinTracks(10)],'FoF',Pole1Ref,'name','PoleRefKinI10Full','crop','full')
project1D(MD,[P1TrackISO,inlierKinTracks(10)],'FoF',Pole1Ref,'name','PoleRefKinI10')
project1D(MD,[P1TrackISO,inlierKinTracks(10)],'FoF',InterpolarRefISO,'name','InterpolarRefKinI10','crop','full');
project1D(MD,[P1TrackISO,inlierKinTracks(10)],'FoF',InterpolarRefISO,'name','InterpolarManifoldCropRefKinI10','crop','manifold');

%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(2)])
%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(5)])
