function [overlayCell]=buildPoleKT3DMovies(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',MD.searchPackageName('buildPoleKT3DMovies'));
    ip.addParameter('trackingKTPackage',[]);
    ip.addParameter('buildSpindleRefPackage',[]);
    ip.addParameter('debug',[]);
    ip.addParameter('name','buildPoleKT3DMovies');
    ip.addParameter('forceRunIdx',[]);
    ip.addParameter('dynROIs',[]);
    ip.addParameter('KTP',[]);
    ip.parse(varargin{:});
    p=ip.Results;

temporaryPackage=GenericPackage({ ... 
    ExternalProcess(MD,'build3DROImovies'),...  
  },[],'name_','buildPoleKT3DMovies-Backup');

% If package with same name exist, retrieve it, otherwise add a new one.
package=[];
packageIdx=0;
if(isempty(p.package))
    package=temporaryPackage;
    MD.addPackage(package);
    packageIdx=length(MD.packages_);
else
    [package]=MD.searchPackageName(p.name);
    packageIdx=package.getIndex();
end

%% Loading some fiber results (not optimal)
if(isempty(p.trackingKTPackage))
    fiberTrackabilityPackage=MD.searchPackageName('fiberTrackabilityAnalysis','selectIdx',1);
    processDetectKT=fiberTrackabilityPackage.getProcess(3);
    processTrackKT=fiberTrackabilityPackage.getProcess(4);
    processDetectPoles=fiberTrackabilityPackage.getProcess(5);
    processBuildRef=fiberTrackabilityPackage.getProcess(6);
else
    processDetectKT=p.trackingKTPackage.getProcess(1);
    processTrackKT=p.trackingKTPackage.getProcess(3);
    processDetectPoles=p.buildSpindleRefPackage.getProcess(1);
    processBuildRef=p.buildSpindleRefPackage.getProcess(2);
end

tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);
fiducials=[P1 P2];

% inlier tracks
refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;

%% Mapping tracks
tmp=load(processTrackKT.outFilePaths_{2}); kinTracksISO=TracksHandle(tmp.tracksFinal);
kinTracksISOInliers=mapTracksTo1DManifold(ROIs{1,2},kinTracksISO,0,'position','start','distType','vertexDistOtsu');

lpid=0;

KTP=p.KTP;
if(isempty(KTP))
    KTP=[(1:100)' ones(100,1);(1:100)' 2*ones(100,1)];
end

disp('Build Lab ref ROIs')
dynROIData=cell(size(KTP,1),1);
for ROIIndex=1:size(KTP,1)
    tIdx=KTP(ROIIndex,1);
    pIdx=KTP(ROIIndex,2);    
    track=kinTracksISOInliers(tIdx);
    PTrack=fiducials(pIdx);
    ROI=[PTrack track];
    ref=FrameOfRef().genCanonicalRef(MD.nFrames_);
    ldynROIData.ref=ref;
    ldynROIData.ROI=ROI;
        % mappedTracks={track.associatedMTP1 track.associatedMTP2};
        % ldynROIData.mappedTracks=mappedTracks;
    dynROIData{ROIIndex}=ldynROIData;
end

% Collect kinROI ROIs
kinROIs=dynROIData;

%% Color KinROIs
% Get Rhos over time and closestPoles;
trackRho=cell(1,length(kinTracksISOInliers));
closestPoleIdx=cell(1,length(kinTracksISOInliers));
for tIdx=1:numel(kinTracksISOInliers)  
    track=kinTracksISOInliers(tIdx);
    trSpindleRefP1P2=refs(1,2).applyBase(track,'');
    trSpindleRefP2P1=refs(2,1).applyBase(track,'');
    [~,~,rhoCellP1P2]=cart2sph(trSpindleRefP1P2.x,trSpindleRefP1P2.y,trSpindleRefP1P2.z);
    [~,~,rhoCellP2P1]=cart2sph(trSpindleRefP2P1.x,trSpindleRefP2P1.y,trSpindleRefP2P1.z);
    [trackRho{tIdx}]=min(rhoCellP2P1,rhoCellP1P2);
    closestPoleIdx{tIdx}=(rhoCellP1P2>rhoCellP2P1)+1;
end

% build1DManifold(trackSet1,trackSet2,varargin)

% Color according to rho
colorIndx=cell(size(kinROIs));
minRho=min(cellfun(@min,trackRho));
maxRho=max(cellfun(@max,trackRho));
for ROIIndex=1:size(KTP,1)
    tIdx=KTP(ROIIndex,1);
    pIdx=KTP(ROIIndex,2);
    colorIndx{ROIIndex}= ceil(127*mat2gray(trackRho{tIdx},[minRho,maxRho])) + 127*(pIdx-1);
end

% Color according to pole
colorIndx=cell(1,numel(kinROIs));
minRho=min(cellfun(@min,trackRho));
maxRho=max(cellfun(@max,trackRho));
for ROIIndex=1:size(KTP,1)
    tIdx=KTP(ROIIndex,1);
    pIdx=KTP(ROIIndex,2);
    colorIndx{ROIIndex}= 127;
end

%% Filter KinROIs
ROIToKeep=ones(1,numel(kinROIs));

% Delete ROI that are linked to the farthest Pole
% endClosestPole=cellfun(@(r) r(end),closestPoleIdx);
% for ROIIndex=1:size(KTP,1)
%     tIdx=KTP(ROIIndex,1);
%     pIdx=KTP(ROIIndex,2);
%     if(pIdx~=endClosestPole(tIdx))  % I write shitty code under stress
%        ROIToDelete(ROIIndex)=1;
%     end
% end

% If KT are very close, delete one (should be 4 ROIs) 
% Keep the one that is the closest to the pole (1 left)
% trackStarts=arrayfun(@(t) [t.x(1) t.y(1) t.z(1)],kinTracksISOInliers,'unif',0);
% trackStarts=vertcat(trackStarts{:});
% D = createSparseDistanceMatrix(trackStarts, trackStarts, 5);
% [sortedScore,indx]=sort(full(D));
% trackStartRhos=cellfun(@(r) r(1),trackRho);
% for ROIIndex=1:size(KTP,1)
%     tIdx=KTP(ROIIndex,1);
%     pIdx=KTP(ROIIndex,2);
%     if(sortedScore(tIdx,end)>0)
%        ROIToDelete(ROIIndex) =(trackStartRhos(tIdx)>trackStartRhos(indx(tIdx,end)));
%     end
% end

%% lifetime
for ROIIndex=1:size(KTP,1)
    tIdx=KTP(ROIIndex,1);
    track=kinTracksISOInliers(tIdx);
    if(track.lifetime<MD.nFrames_)
        ROIToKeep(ROIIndex)=0;
    end 
end

%% Keep the first two ROI because they actually look pretty nice
keepIdx=find(ROIToKeep);
ROIToKeep=zeros(size(ROIToKeep));
ROIToKeep(keepIdx(1:2))=1;

%% Keep only 4 filtered ROI, make sure that the KT are not the same.
% pIdx=KTP(:,2);
% keptROI=(find(ROIToDelete==0));
% % randKept=keptROI(randperm(length(keptROI)));
% keepROIP1=keptROI(pIdx(keptROI)==1);
% keepROIP2=keptROI(pIdx(keptROI)==2);
% ROIToDelete=ones(size(ROIToDelete));
% ROIToDelete([keepROIP1(1) keepROIP2]=0;

%%
kinROIs(~logical(ROIToKeep))=[];
colorIndx(~logical(ROIToKeep))=[];

if(~isempty(p.dynROIs))
        kinROIs=arrayfun(@(r) r,p.dynROIs,'unif',0);
        colorIndx=arrayfun(@(r) 120,p.dynROIs,'unif',0);

end

try
    MD.deletePackage(MD.searchPackageName('ROIVolumes'));
catch
end;
MDCropROI=build3DROImovies(MD,kinROIs,colorIndx,p.name);
printMIP(MDCropROI);
processRenderer1=ExternalProcess(MD,'renderer');
processProj=ExternalProcess(MD,'dynROIProj');
projectDynROI(MD,'processRenderer',processRenderer1,'processSingleProj',processProj);
processRenderer2=ExternalProcess(MDCropROI,'renderer');
projectDynROI(MDCropROI,'processRenderer',processRenderer2,'processSingleProj',processProj);

printProcMIPArray({processRenderer1,processRenderer2},[MDCropROI.outputDirectory_ filesep 'ColorAndOrig'])


completePackage=GenericPackage(temporaryPackage.processes_,[],'name_',p.name);
MD.setPackage(packageIdx,completePackage);

