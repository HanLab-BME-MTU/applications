function [kinTracksInliers,randKinTracksInliers]=MTTipsBias(MD,varargin)
% Provided a spindle pole referential (process 'poleRef') compare the bias
% of +TIPs detection toward the kinetocore referential or a random location
% inside the spindle. 

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addOptional('kinTracksWithSpindle',[]);
ip.addOptional('EB3TracksWithSpindle',[]);
ip.addOptional('randKinTracksWithSpindle',[]);
ip.addOptional('detLabRef',[]);
ip.addParameter('kinRange',[],@isnumeric);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('cutoffs',250, @isnumeric);
ip.addParameter('distType','angle');
ip.addParameter('plotHandle',[]);
ip.addParameter('process',[]);
ip.parse(MD,varargin{:});
p=ip.Results;
cutoff=p.cutoffs;
kinRange=p.kinRange;

printAll=p.printAll;

%% Loading kin and EB3
poleRefProcessIdx=(cellfun(@(p) strcmp(p.name_,'poleRef'),MD.processes_));
if(any(poleRefProcessIdx)&&isempty(p.kinTracksWithSpindle))
    disp('reading process');tic;
    poleRefProcess=MD.getProcess(find(poleRefProcessIdx,1,'last'));
    EB3Tracks=load(poleRefProcess.outFilePaths_{1}); EB3Tracks=EB3Tracks.EB3Tracks;
    kinTracks=load(poleRefProcess.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
    EB3TracksInliers=load(poleRefProcess.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
    kinTracksInliers=load(poleRefProcess.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;
    detLabRef=load(poleRefProcess.outFilePaths_{5}); detLabRef=detLabRef.detLabRef;
    toc
else
    EB3Tracks=p.EB3TracksWithSpindle;
    kinTracks=p.kinTracksWithSpindle;
    inliersEB3=(logical(arrayfun(@(eb) eb.inliers(1),EB3Tracks)));
    EB3TracksInliers=EB3Tracks(inliersEB3);
    inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracks));
    kinTracksInliers=kinTracks(inliersKin);
    detLabRef=p.detLabRef;
end

%% Loading rand kin 
if(isempty(kinRange))
    kinRange=1:length(kinTracksInliers)
end

%% Functionalize and processize
randKinProcIdx=cellfun(@(p) strcmp(p.name_,'randKin'),MD.processes_);
if(any(randKinProcIdx)&&isempty(p.randKinTracksWithSpindle))
    randKinProc=MD.getProcess(find(randKinProcIdx,1,'last'));
    tmp=load(randKinProc.outFilePaths_{1}); 
    randKinTracks=tmp.randKinTracks;
else
    randKinTracks=p.randKinTracksWithSpindle;
end
%%

inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracks));
randKinTracksInliers=randKinTracks(inliersKin);

dummyProc=ExternalProcess(MD);
% Bias with detection
% mapDetectionsToKin(kinTracksInliers(kinRange),detLabRef,cutoff)
% appearingMTRefProjection(kinTracksInliers(kinRange),'associatedDetectP1','associatedDetectP2');
% printAllDetectKinPoleRef(kinTracksInliers(kinRange),'name','detectBased','process',dummyProc)
% 
% mapDetectionsToKin(randKinTracksInliers(kinRange),detLabRef,cutoff)
% appearingMTRefProjection(randKinTracksInliers(kinRange),'associatedDetectP1','associatedDetectP2');
% printAllDetectKinPoleRef(randKinTracksInliers(kinRange),'name','randDetectBased','process',dummyProc)

% Bias with tracked detection
mapMTTipsToKin(kinTracksInliers(kinRange),EB3TracksInliers,cutoff)
appearingMTRefProjection(kinTracksInliers(kinRange),'associatedTipsP1','associatedTipsP2');
printAllMTTipsKinPoleRef(kinTracksInliers(kinRange),'name','trackBased','process',dummyProc)

mapMTTipsToKin(randKinTracksInliers(kinRange),EB3TracksInliers,cutoff)
appearingMTRefProjection(randKinTracksInliers(kinRange),'associatedTipsP1','associatedTipsP2');
printAllMTTipsKinPoleRef(randKinTracksInliers(kinRange),'name','randTrackBased','process',dummyProc)

MTTipsBiasProcess=p.process;
%%
procFolder=[MD.outputDirectory_  filesep 'Kin' filesep 'MTTipsBias' filesep];
mkdir(procFolder);
if(~isempty(MTTipsBiasProcess))
    save([procFolder 'MTTipsBiasForKinAndRand.mat'],'kinTracksInliers','randKinTracksInliers')
    MTTipsBiasProcess.setOutFilePaths({[procFolder 'MTTipsBiasForKinAndRand.mat']})
end;
%%
kIdx=1;ktInlier=kinTracksInliers(kIdx);rktInlier=randKinTracksInliers(kIdx);

[Hs]=setupFigure(1,2,2,'AxesWidth',8,'AxesHeight',4);
printMTTipsPoleRef(ktInlier.KP1,ktInlier.associatedTipsP2KinRef,ktInlier.associatedTipsP2Idx,'handle',Hs(1))
printMTTipsPoleRef(rktInlier.KP1,rktInlier.associatedTipsP2KinRef,rktInlier.associatedTipsP2Idx,'handle',Hs(2))

%%
displayTipsBiasStat({kinTracksInliers(kinRange),randKinTracksInliers(kinRange)},{'kin','rand'});

disp('end');
