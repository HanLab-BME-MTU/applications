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
P1Test=TracksHandle();
P1Test.endFrame=20;
P1Test.startFrame=1;
P1Test.x=arrayfun(@(d) d.xCoord(1,1),poleMovieInfo)';
P1Test.y=arrayfun(@(d) d.yCoord(1,1),poleMovieInfo)';
P1Test.z=arrayfun(@(d) d.zCoord(1,1)*MD.pixelSizeZ_/MD.pixelSize_,poleMovieInfo)'

project1D(MD,[P1Test,inlierKinTracks(1)])
%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(2)])
%project1D(MD,[P1Test,kinTracksRawData.tracksLabRef(5)])