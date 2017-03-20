function randKinTracks=randomizeKinSpindleRef(MD,randomDist,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('process',[]);
ip.parse(MD,varargin{:});
p=ip.Results;
% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
kinTracksOrig=kinTrackData.tracksLabRef;
[randKinTracks]=randomizeKinetochore(kinTracksOrig,randomDist);

% Translate these changes in the detection structure and associated polar
% coordiante
% Load associated data
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionsLabRef=tmp.detectionsLabRef;

dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
EB3poleDetectionMethod=['simplex_scale_003'];
outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep EB3poleDetectionMethod filesep];
poleData=load([outputDirPoleDetect filesep 'poleDetection.mat']);
poleMovieInfo=poleData.poleMovieInfo;

[~,kinSphericalCoord,inliers]=tracks2detections(randKinTracks,detectionsLabRef,poleMovieInfo,dataIsotropy);

% Rebuild the augmented kin
[randKinTracks]=addSpindleRef(MD,'kinTracks',randKinTracks,'kinSphericalCoord',kinSphericalCoord,'kinInliers',inliers);

process=p.process;
if(~isempty(process))
    %%
    procFolder=[MD.outputDirectory_  filesep 'Kin' filesep 'randomized' filesep];
    mkdir(procFolder);
    save([procFolder 'randKinTracks.mat'],'randKinTracks')
    process.setOutFilePaths({[procFolder 'randKinTracks.mat']})
end;

end
