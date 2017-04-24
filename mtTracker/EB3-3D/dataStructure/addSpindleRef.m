function [kinTracks,EB3Tracks,detLabRef]=addSpindleRef(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('MD',[],@(MD) isa(MD,'MovieData'));
ip.addParameter('kinTracks',[],@(x) isa(x,'Tracks'));
ip.addParameter('kinSphericalCoord',[]);
ip.addParameter('kinInliers',[]);
ip.addParameter('processDetectPoles',[]);
ip.addParameter('EB3tracks',[],@(x) isa(x,'Tracks'));
ip.addParameter('EB3SphCoord',[]);
ip.addParameter('EB3Inlier',[]);
ip.addParameter('EB3PoleId',[]);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('process',[]);
ip.addParameter('testKinIdx',[19 46 156],@isnumeric);
ip.parse(MD,varargin{:});
p=ip.Results;

% Augment the structures with spherical Coordinate.
%% Load the pole info
MD=p.MD;
if(~isempty(p.processDetectPoles))
    tmp=load(p.processDetectPoles.outFilePaths_{1});
    poleMovieInfo=tmp.poleMovieInfo;
    MD=p.processDetectPoles.getOwner();
end


poleRefs=buildPoleRef(poleMovieInfo,MD.pixelSize_);

%% Load EB3 tracks add azimuth info, change coordinate to real space measurement.
if(isempty(p.EB3tracks)||isempty(p.EB3SphCoord))
    outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'detection' filesep];
    tmp=load([outputDirDetect filesep 'sphericalCoord.mat']);
    EB3SphCoord=tmp.sphCoord;
    tmp=load([outputDirDetect filesep 'detectionLabRef.mat']);
    EB3LabRef=tmp.detectionsLabRef;
    EB3PoleDist=load([outputDirDetect filesep 'dist.mat']);
    EB3PoleId=EB3PoleDist.poleId;
    EB3Inliers=EB3PoleDist.inliers;
    outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
    tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
    EB3Tracks=tmp.tracksLabRef;

    dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
    [~,EB3SphCoord,EB3PoleId,EB3Inliers,~,~]=poleDist(poleMovieInfo,EB3LabRef,'anisotropy',dataIsotropy,'angleRef','poles');
else
    EB3Tracks=p.EB3tracks;
    EB3SphCoord=p.EB3SphCoord;
    EB3Inliers=p.EB3Inliers;
    EB3PoleId=p.EB3PoleId;
end

EB3Tracks=addSpindleRefEB3(MD,poleRefs,EB3Tracks,EB3SphCoord,EB3Inliers,EB3PoleId);

%%
if(isempty(p.kinTracks)||isempty(p.kinSphericalCoord))
    outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
    kinSphericalCoord=load([outputDirDetect filesep 'sphericalCoordBothPoles.mat']);
    kinSphericalCoord=kinSphericalCoord.sphCoord;
    kinPoleDist=load([outputDirDetect filesep 'dist.mat']);
    kinInliers=kinPoleDist.inliers;
    outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
    kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
    kinTracks=kinTrackData.tracksLabRef;

    outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];mkdir(outputDirDetect);
    detectionsLabRef=load([outputDirDetect filesep 'detectionLabRef.mat']);
    detectionsLabRef=detectionsLabRef.detectionsLabRef;

    dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
    [~,kinSphericalCoord,~,kinInliers,~,~]=poleDist(poleMovieInfo,detectionsLabRef,'anisotropy',dataIsotropy,'angleRef','poles');
else
    kinTracks=p.kinTracks;
    kinSphericalCoord=p.kinSphericalCoord;
    kinInliers=p.kinInliers;
end

kinTracks=addSpindleRefKin(MD,poleRefs,kinTracks,kinSphericalCoord,kinInliers);
% If there is a process object, write data directly in appropriate folder
if(~isempty(p.process))

    poleRefProcess=p.process;
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'Kin' filesep 'track'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'kinTracks');
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'EB3' filesep 'track'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'EB3Tracks');
%     outputDirCatchingMT=[MD.outputDirectory_ filesep 'EB3' filesep 'detection'];
%     save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'detLabRef');
    save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'EB3Tracks');


    %% inlier index
    inliersEB3=(logical(arrayfun(@(eb) eb.inliers(1),EB3Tracks)));
    inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracks));
    kinTracksInliers=kinTracks(inliersKin);
    EB3TracksInliers=EB3Tracks(inliersEB3);

    outputDirCatchingMT=[MD.outputDirectory_ filesep 'Kin' filesep 'track'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRefInliers.mat'],'kinTracksInliers');
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'EB3' filesep 'track'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRefInliers.mat'],'EB3TracksInliers');
    %%
    poleRefProcess.setOutFilePaths({ [MD.outputDirectory_ '' filesep 'EB3' filesep 'track' filesep 'augmentedSpindleRef.mat'], ...
        [MD.outputDirectory_ '' filesep 'Kin' filesep 'track' filesep 'augmentedSpindleRef.mat'], ...
        [MD.outputDirectory_ '' filesep 'EB3' filesep 'track' filesep 'augmentedSpindleRefInliers.mat'], ...
        [MD.outputDirectory_ '' filesep 'Kin' filesep 'track' filesep 'augmentedSpindleRefInliers.mat'],  ...

%         [MD.outputDirectory_ '' filesep 'EB3' filesep 'detection' filesep 'augmentedSpindleRef.mat'] ...
        });
    pa = poleRefProcess.getParameters();
    pa.parameters = ip.Results;
    poleRefProcess.setParameters(pa);
    poleRefProcess.setDateTime();

end
%%
