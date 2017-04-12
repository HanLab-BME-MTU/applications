function randKinTracks=randomizeKinSpindleRef(MD,randomDist,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('randomDist');
ip.addParameter('process',[]);
ip.addParameter('processDetectPoles',[]);

ip.parse(MD,randomDist,varargin{:});
p=ip.Results;

% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
kinTracksOrig=kinTrackData.tracksLabRef;
[randKinTracks]=randomizeKinetochore(kinTracksOrig,randomDist);

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

dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
[~,kinSphericalCoord,inliers]=tracks2detections(randKinTracks,detectionsLabRef,poleMovieInfo,dataIsotropy);

% Rebuild the augmented kin
randKinTracks=addSpindleRefKin(MD,poleMovieInfo,randKinTracks,kinSphericalCoord,inliers);

process=p.process;
if(~isempty(process))
    %%
    procFolder=[MD.outputDirectory_  filesep 'Kin' filesep 'randomized' filesep];
    mkdir(procFolder);
    save([procFolder 'randKinTracks.mat'],'randKinTracks');
    process.setOutFilePaths({[procFolder 'randKinTracks.mat']})
end;

end
