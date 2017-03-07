function [kinTracksInliers,randKinTracksInlier]=MTTipsBias(MD,varargin)
% Provided a spindle pole referential (process 'poleRef') compare the bias
% of +TIPs detection toward the kinetocore referential or a random location
% inside the spindle. 

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('kinTracksWithSpindle',[]);
ip.addParameter('randKinTracksWithSpindle',[]);
ip.addParameter('EB3TracksWithSpindle',[]);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('cutoffs',250, @isnumeric);
ip.addParameter('distType','angle');
ip.addParameter('plotHandle',[]);
ip.addParameter('process',[]);
ip.addParameter('rerandomize',false);
ip.parse(MD,varargin{:});
p=ip.Results;

printAll=p.printAll;
%%
EB3Tracks=p.EB3TracksWithSpindle;
kinTracks=p.kinTracksWithSpindle;
poleRefProcessIdx=(cellfun(@(p) strcmp(p.name_,'poleRef'),MD.processes_));
if(any(poleRefProcessIdx))
    poleRefProcess=MD.getProcess(find(poleRefProcessIdx,1,'last'));
    EB3Tracks=load(poleRefProcess.outFilePaths_{1}); EB3Tracks=EB3Tracks.EB3Tracks;
    kinTracks=load(poleRefProcess.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
    EB3TracksInliers=load(poleRefProcess.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
    kinTracksInliers=load(poleRefProcess.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;
    detLabRef=load(poleRefProcess.outFilePaths_{5}); detLabRef=detLabRef.detLabRef;
end
%%


% TODO:  functionalize and processize
randKinProcIdx=cellfun(@(p) strcmp(p.name_,'randKin'),MD.processes_);
if(any(randKinProcIdx)&&~(p.rerandomize))
    randKinProc=MD.getProcess(find(randKinProcIdx,1,'last'));
    tmp=load(randKinProc.outFilePaths_{1}); 
    randKinTracks=tmp.randKinTracks;
    randKinTracksInlier=tmp.randKinTracksInlier;
else
    if(isempty(p.randKinTracksWithSpindle))
        %%
        % Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
        randomDist=10; % in pixel
        
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
        [randKinTracksPlus]=addSpindleRef(MD,'kinTracks',randKinTracks,'kinSphericalCoord',kinSphericalCoord,'kinInliers',inliers);
        
        % Estimate bundle outside the Kin-Plan refencial
        inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracks));
        randKinTracksInlier=randKinTracksPlus(inliersKin);
        %%
    else
        randKinTracksInlier=p.randKinTracksWithSpindle;
    end
end

cutoff=p.cutoffs;
dummyProc=ExternalProcess(MD);
% Bias with detection
mapDetectionsToKin(kinTracksInliers,detLabRef,cutoff)
appearingMTRefProjection(kinTracksInliers,'associatedDetectP1','associatedDetectP2');
printAllDetectKinPoleRef(kinTracksInliers,'name','detectBased','process',dummyProc)

mapDetectionsToKin(randKinTracksInlier,detLabRef,cutoff)
appearingMTRefProjection(randKinTracksInlier,'associatedDetectP1','associatedDetectP2');
printAllDetectKinPoleRef(randKinTracksInlier,'name','detectBased','process',dummyProc)

% Bias with tracked detection
mapMTTipsToKin(kinTracksInliers,EB3TracksInliers,cutoff)
appearingMTRefProjection(kinTracksInliers,'associatedTipsP1','associatedTipsP2');
printAllMTTipsKinPoleRef(kinTracksInliers,'name','trackBased','process',dummyProc)

mapMTTipsToKin(randKinTracksInlier,EB3TracksInliers,cutoff)
appearingMTRefProjection(randKinTracksInlier,'associatedTipsP1','associatedTipsP2');
printAllMTTipsKinPoleRef(randKinTracksInlier,'name','randTrackBased','process',dummyProc)

MTTipsBiasProcess=p.process;
%%
procFolder=[MD.outputDirectory_ 'Kin' filesep 'MTTipsBias'];
mkdir(procFolder);
if(~isempty(MTTipsBiasProcess))
    save([procFolder 'MTTipsBiasForKinAndRand.mat'],'kinTracksInliers','randKinTracksInlier')
    MTTipsBiasProcess.setOutFilePaths({[procFolder 'MTTipsBiasForKinAndRand.mat']})
end;
%%
kIdx=1;ktInlier=kinTracksInliers(kIdx);rktInlier=randKinTracksInlier(kIdx);

[Hs]=setupFigure(1,2,2,'AxesWidth',8,'AxesHeight',4);
printMTTipsPoleRef(ktInlier.KP1,ktInlier.associatedTipsP2KinRef,ktInlier.associatedTipsP2Idx,'handle',Hs(1))
printMTTipsPoleRef(rktInlier.KP1,rktInlier.associatedTipsP2KinRef,rktInlier.associatedTipsP2Idx,'handle',Hs(2))

%%
displayTipsBiasStat({kinTracksInliers,randKinTracksInlier},{'kin','rand'});

disp('end');
