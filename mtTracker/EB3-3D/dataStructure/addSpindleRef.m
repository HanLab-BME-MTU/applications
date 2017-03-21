function [kinTracks,EB3Tracks,detLabRef]=addSpindleRef(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('MD',[],@(MD) isa(MD,'MovieData'));
ip.addParameter('kinTracks',[],@(x) isa(x,'Tracks'));
ip.addParameter('processDetectPoles',[]);
ip.addParameter('EB3tracks',[],@(x) isa(x,'Tracks'));
ip.addParameter('kinSphericalCoord',[]);
ip.addParameter('EB3SphCoord',[]);
ip.addParameter('EB3Inlier',[]);
ip.addParameter('EB3PoleId',[]);
ip.addParameter('kinInliers',[]);
ip.addParameter('kinPoleId',[]);
ip.addParameter('kinPoleDist',[]);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('process',[]);
ip.addParameter('processInput',[]);
ip.addParameter('testKinIdx',[19 46 156],@isnumeric);
ip.addParameter('distanceCutOff',0.1,@isnumeric);
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

% WARNING: this is not a trajectory, merely a collection of poles to ease
% implementation.
P1=struct();
P1.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(1,1)-1)+1,poleMovieInfo)';
P1.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(1,1)-1)+1,poleMovieInfo)';
P1.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(1,1)-1)+1,poleMovieInfo)';

P2=struct();
P2.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(2,1)-1)+1,poleMovieInfo)';
P2.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(2,1)-1)+1,poleMovieInfo)';
P2.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(2,1)-1)+1,poleMovieInfo)';

refP1=FrameOfRef();
refP1.setOriginFromTrack(P1);
refP1.setZFromTrack(P2);
refP1.genBaseFromZ();

refP2=FrameOfRef();
refP2.setOriginFromTrack(P2);
refP2.setZFromTrack(P1);
refP2.genBaseFromZ();

poleRefs=[refP1 refP2];

% For MT detection
tic;
detLabRef=Detections(EB3LabRef);
detLabRef.scale(MD);
dp1=refP1.applyBaseToDetection(detLabRef,'pole1');
dp2=refP2.applyBaseToDetection(detLabRef,'pole2');
dp1.addSphericalCoord();
dp2.addSphericalCoord();
toc

%% For MT tracks
tic;
for tIdx=1:length(EB3Tracks)
    %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');

    tr=EB3Tracks(tIdx);
     try
        tr.addprop('inliers');
        tr.addprop('poleId');
        tr.addprop('azimuth');      % DEPRECATED
        tr.addprop('elevation');    % DEPRECATED
        tr.addprop('rho');          % DEPRECATED      
    catch
    end;
    tr.x=(tr.x-1)*MD.pixelSize_+1;
    tr.y=(tr.y-1)*MD.pixelSize_+1;
    tr.z=(tr.z-1)*MD.pixelSize_+1;

    nonGap=~tr.gapMask();
    tr.poleId=nan(size(tr.f));
    tr.poleId(nonGap)=arrayfun(@(i,f) EB3PoleId{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.inliers=nan(size(tr.f));
    tr.inliers(nonGap)=arrayfun(@(i,f) EB3Inliers{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));    
    
    % DEPRECATED
    tr.azimuth=nan(2,length(tr.f));
    tr.elevation=nan(2,length(tr.f));
    tr.rho=nan(2,length(tr.f));
    for poleID=1:2
        tr.azimuth(poleID,nonGap)=arrayfun(@(i,f) EB3SphCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.elevation(poleID,nonGap)=arrayfun(@(i,f) EB3SphCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.rho(poleID,nonGap)=arrayfun(@(i,f) EB3SphCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    end 
    % END DEPRECATED
end
toc;
%%
tic;
refP1.applyBaseToTrack(EB3Tracks,'pole1');
refP2.applyBaseToTrack(EB3Tracks,'pole2');
%%
% augment pole ref projected data with spherical coordinate
refName={'pole1','pole2'};
for tIdx=1:length(EB3Tracks)
%     trackPoleRefs=[];
    for poleID=1:length(poleRefs);
        % Copying EB3 track
        tr=getfield(EB3Tracks(tIdx),refName{poleID});
        % Adding correponding spherical coordinate
        try
            tr.addprop('azimuth');
            tr.addprop('elevation');
            tr.addprop('rho');
        catch
        end;
               
        nonGap=~tr.gapMask();
        tr.azimuth=nan(1,length(tr.f));
        tr.elevation=nan(1,length(tr.f));
        tr.rho=nan(1,length(tr.f));       
         
        tr.azimuth(nonGap)=arrayfun(@(i,f) EB3SphCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.elevation(nonGap)=arrayfun(@(i,f) EB3SphCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.rho(nonGap)=arrayfun(@(i,f) EB3SphCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    end
end
toc

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

tic;
% Augment the structures with spherical Coordinate. 
for kIdx=1:length(kinTracks)
    %progressText(kIdx/length(kinTracks),'Loading kin spherical coordinates.');
    tr=kinTracks(kIdx);
    try
        tr.addprop('inliers');
        tr.addprop('azimuth');
        tr.addprop('elevation');
        tr.addprop('rho');
    catch
    end;
    tr.x=(tr.x-1)*MD.pixelSize_+1;
    tr.y=(tr.y-1)*MD.pixelSize_+1;
    tr.z=(tr.z-1)*MD.pixelSize_+1;
    
    nonGap=~tr.gapMask();
    tr.azimuth=nan(2,length(tr.f));
    tr.elevation=nan(2,length(tr.f));
    tr.rho=nan(2,length(tr.f));
    tr.inliers=nan(size(tr.f));
    tr.inliers(nonGap)=arrayfun(@(i,f) kinInliers{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    for poleID=1:2
        tr.azimuth(poleID,nonGap)=arrayfun(@(i,f) kinSphericalCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.elevation(poleID,nonGap)=arrayfun(@(i,f) kinSphericalCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.rho(poleID,nonGap)=arrayfun(@(i,f) kinSphericalCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    end 
end
%%
toc;
tic;
refP1.applyBaseToTrack(kinTracks,'pole1');
refP2.applyBaseToTrack(kinTracks,'pole2');

for tIdx=1:length(kinTracks)
%     trackPoleRefs=[];
    for poleID=1:length(poleRefs);
        % Copying EB3 track
        tr=getfield(kinTracks(tIdx),refName{poleID});
        % Adding correponding spherical coordinate
        try
            tr.addprop('azimuth');
            tr.addprop('elevation');
            tr.addprop('rho');
        catch
        end;
               
        nonGap=~tr.gapMask();
        tr.azimuth=nan(1,length(tr.f));
        tr.elevation=nan(1,length(tr.f));
        tr.rho=nan(1,length(tr.f));       
         
        tr.azimuth(nonGap)=arrayfun(@(i,f) kinSphericalCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.elevation(nonGap)=arrayfun(@(i,f) kinSphericalCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.rho(nonGap)=arrayfun(@(i,f) kinSphericalCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    end
end
toc

% If there is a process object, write data directly in appropriate folder
if(~isempty(p.process))

    poleRefProcess=p.process;
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'Kin' filesep 'track'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'kinTracks');
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'EB3' filesep 'track'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'EB3Tracks');
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'EB3' filesep 'detection'];
    save([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'detLabRef');
    
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
        [MD.outputDirectory_ '' filesep 'EB3' filesep 'detection' filesep 'augmentedSpindleRef.mat'] ...
        });
    pa = poleRefProcess.getParameters();
    pa.parameters = ip.Results;
    poleRefProcess.setParameters(pa);
    poleRefProcess.setDateTime();

end
%%

