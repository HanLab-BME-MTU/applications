function [randKinTracks,randKinTracksISO]=randomizeDenseTrackBruteForce(track,manifolds,maxRandomDist,varargin)
% randomize a track so that no tracks map in the manifolds
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('maxRandomDist');
ip.addRequired('mappingMetric');
ip.addParameter('process',[]);
ip.addParameter('processDetectPoles',[]);
ip.parse(MD,maxRandomDist,varargin{:});
p=ip.Results;

% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
kinTracksOrig=kinTrackData.tracksLabRef;
[randKinTracksISO]=randomizeKinetochore(kinTracksOrig,maxRandomDist);

% saturating randomization with the volume limits
for tIdx=1:length(randKinTracksISO)
  randKin=randKinTracksISO(tIdx)
  randKin.x=min(randKin.x,MD.imSize_(1));
  randKin.x=max(randKin.x,1);
  randKin.y=min(randKin.y,MD.imSize_(2));
  randKin.y=max(randKin.y,1);
  randKin.z=min(randKin.z,MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_);
  randKin.z=max(randKin.z,1);
end
% Translate these changes in the detection structure and associated polar
% coordiante
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionsLabRef=tmp.detectionsLabRef;


% Load associated data
if(~isempty(p.processDetectPoles))
    tmp=load(p.processDetectPoles.outFilePaths_{1});
    poleMovieInfo=tmp.poleMovieInfo;
end

% Soon to be deprecated (when the code will be fully manifold orientied)
[detections]=tracks2detections(randKinTracksISO,detectionsLabRef);
dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
[~,kinSphericalCoord,~,inliers]=poleDist(poleMovieInfo,detections,'anisotropy',dataIsotropy,'angleRef','poles');

% Build the associated polar coordinate
% WARNING: this is not a trajectory, merely a collection of poles to ease
% implementation.
P1=struct();
P1.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(1,1)-1)+1,poleMovieInfo)';
P1.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(1,1)-1)+1,poleMovieInfo)';
P1.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(1,1)-1)+1,poleMovieInfo)';
P1.f=1:length(poleMovieInfo);

P2=struct();
P2.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(2,1)-1)+1,poleMovieInfo)';
P2.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(2,1)-1)+1,poleMovieInfo)';
P2.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(2,1)-1)+1,poleMovieInfo)';
P2.f=1:length(poleMovieInfo);


refP1=FrameOfRef();
refP1.setOriginFromTrack(P1);
refP1.setZFromTrack(P2);
refP1.genBaseFromZ();

refP2=FrameOfRef();
refP2.setOriginFromTrack(P2);
refP2.setZFromTrack(P1);
refP2.genBaseFromZ();

poleRefs=[refP1 refP2];

% rebuild the augmented kin strucure
randKinTracks=randKinTracksISO.copy();
randKinTracks=addSpindleRefKin(MD,poleRefs,randKinTracks,kinSphericalCoord,inliers);

process=p.process;
if(~isempty(process))
    %%
    procFolder=[MD.outputDirectory_  filesep 'Kin' filesep 'randomized' filesep];
    mkdir(procFolder);
    save([procFolder 'randKinTracks.mat'],'randKinTracks','randKinTracksISO');
    process.setOutFilePaths({[procFolder 'randKinTracks.mat']})
    pa = process.getParameters();
    pa = ip.Results;
    process.setParameters(pa);
    process.setDateTime();
end;

end
