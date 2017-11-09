function bundleStatisticsWholeSpindleTestingInlier(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('printAll',false, @islogical);
ip.addParameter('testKinIdx',[19 46 156],@isnumeric);
ip.parse(MD,varargin{:});
p=ip.Results;

outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];


[kinTracksPlus,EB3tracksPlus]=addSpindleRef(MD);

%% Estimate bundle in KinPole axis
captureDetection(MD,'name','allTracks');
allFiber=detectMTAlignment(MD,'name','allTracks');

EB3tracksInliers=EB3tracksPlus(logical(arrayfun(@(eb) eb.inliers(1),EB3tracksPlus)));
inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracksPlus));
kinTracksInliers=kinTracksPlus(inliersKin);

captureDetection(MD,'kinTracks',kinTracksInliers,'EB3tracks',EB3tracksInliers,'name','inliers');
allFiberInlier=detectMTAlignment(MD,'name','inliers');

%% For each kinetochore, plot an Amira file with attached mt
outputDirAmira=[outputDirBundle filesep 'Amira_normal' filesep];
parfor kIdx=1:length(allFiberInlier)
    kinTrack=allFiberInlier(kIdx);
    trackSet=[kinTrack; kinTrack.catchingMT];
    trackType=[1; zeros(length(kinTrack.catchingMT),1)];
    bundleInfo=[0 kinTrack.fiber+1];
    amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType},{'bundle',bundleInfo}})
end
%%
randomDist=10; % in pixel

% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
kinTracks=kinTrackData.tracksLabRef;
[randKinTracks]=randomizeKinetochore(kinTracks,randomDist);

% Translate these changes in the detection structure and associated polar
% coordiante
% Load associated data
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep]
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionsLabRef=tmp.detectionsLabRef;

dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
EB3poleDetectionMethod=['simplex_scale_003'];
outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep EB3poleDetectionMethod filesep];
poleData=load([outputDirPoleDetect filesep 'poleDetection.mat']);
poleMovieInfo=poleData.poleMovieInfo;

[~,kinSphericalCoord,inliers]=tracks2detections(randKinTracks,detectionsLabRef,poleMovieInfo,dataIsotropy)

% Rebuild the augmented kin
[randKinTracksPlus]=addSpindleRef(MD,'kinTracks',randKinTracks,'kinSphericalCoord',kinSphericalCoord,'kinInliers',inliers);

if(p.printAll)
    amiraWriteTracks([MD.outputDirectory_ filesep 'Kin' filesep 'randomized' filesep 'random.am'], randKinTracks ,'scales',dataIsotropy)
    amiraWriteTracks([MD.outputDirectory_ filesep 'Kin' filesep 'randomizedAugmented' filesep 'random.am'], randKinTracksPlus)
end

% % Estimate bundle outside the Kin-Plan refencial
captureDetection(MD,'kinTracks',randKinTracksPlus,'EB3tracks',EB3tracksPlus,'name','allTracksRandom');
allFiberRandom=detectMTAlignment(MD,'name','allTracksRandom');

randKinTracksInlier=randKinTracks(inliersKin);
captureDetection(MD,'kinTracks',randKinTracksInlier,'EB3tracks',EB3tracksInliers,'name','inliersRandom');
allFiberRandomInlier=detectMTAlignment(MD,'name','inliersRandom');
%% For each kinetochore, plot an Amira file with attached mt
outputDirAmira=[outputDirBundle filesep 'Amira_random' filesep];
parfor kIdx=1:length(allFiberRandomInlier)
    kinTrack=allFiberRandomInlier(kIdx);
    trackSet=[kinTrack; kinTrack.catchingMT];
    trackType=[1; zeros(length(kinTrack.catchingMT),1)];
    bundleInfo=[0 kinTrack.fiber+1];
    amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType},{'bundle',bundleInfo}})
end

bundleStatistics(MD,'kinBundle',{allFiber,allFiberInlier,allFiberRandom,allFiberRandomInlier},'kinBundleName',{'allFiber','allFiberInlier','allFiberRandom','allFiberRandomInlier'});



