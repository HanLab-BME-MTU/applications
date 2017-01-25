function bundleStatisticsWholeSpindle(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('printAll',false, @islogical);
ip.addParameter('testKinIdx',[19 46 156],@isnumeric);
ip.parse(MD,varargin{:});
p=ip.Results;

randomDist=10;

% Estimate bundle in KinPole axis
captureDetection(MD);
detectMTAlignment(MD);
bundleStatistics(MD);

% Keep inliers only
outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'detection' filesep];
tmp=load([outputDirDetect filesep 'sphericalCoord.mat']);
EB3SphCoord=tmp.sphCoord;
EB3PoleDist=load([outputDirDetect filesep 'dist.mat']);

outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EB3tracks=tmp.tracksLabRef;

% Randomize Kinetochore and create the associted sphercial coordinates.
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
kinTracks=kinTrackData.tracksLabRef;
[randKinTracks]=randomizeKinetochore(kinTracks,randomDist);

% Translate these changes in the detection structure
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionsLabRef=tmp.detectionsLabRef;
for kIdx=1:length(randKinTracks)
    tr=randKinTracks(kIdx);
    for tIdx=1:tr.lifetime
        f=tr.f(tIdx);
        if(tr.tracksFeatIndxCG(tIdx)>0)
            detectionsLabRef(f).xCoord(tr.tracksFeatIndxCG(tIdx),1)=tr.x(tIdx);
            detectionsLabRef(f).yCoord(tr.tracksFeatIndxCG(tIdx),1)=tr.y(tIdx);
            detectionsLabRef(f).zCoord(tr.tracksFeatIndxCG(tIdx),1)=tr.z(tIdx);
        end
   end
end

% Build the associated polar coordinate
dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
EB3poleDetectionMethod=['simplex_scale_003'];
outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep EB3poleDetectionMethod filesep];
poleData=load([outputDirPoleDetect filesep 'poleDetection.mat']);
poleMovieInfo=poleData.poleMovieInfo;
[~,kinSphericalCoord]=poleDist(poleMovieInfo,detectionsLabRef,'anisotropy',dataIsotropy,'angleRef','poles');

%[randKinTracks]=randomizeKinetochore(kinTracks,randomDist);
if(p.printAll)
amiraWriteTracks([MD.outputDirectory_ filesep 'Kin' filesep 'randomized' filesep 'random.am'], randKinTracks ,'scales',dataIsotropy)
end

% Estimate bundle outside the Kin-Plan refencial
captureDetection(MD,'kinTracks',randKinTracks,'kinSphericalCoord',kinSphericalCoord)
detectMTAlignment(MD)
bundleStatistics(MD)