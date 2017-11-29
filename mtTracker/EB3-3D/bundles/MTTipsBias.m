function [kinTracksInliers,randKinTracksInliers]=MTTipsBias(varargin)
% Provided a spindle pole referential (process 'poleRef') compare the bias
% of +TIPs detection toward the kinetocore referential or a random location
% inside the spindle. 

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('processSpindleRef',[]);
ip.addOptional('processRandKin',[]);
ip.addOptional('kinTracksWithSpindle',[]);
ip.addOptional('EB3TracksWithSpindle',[]);
ip.addOptional('randKinTracksWithSpindle',[]);
ip.addOptional('detLabRef',[]);
ip.addParameter('process',[]);
ip.addParameter('kinRange',[],@isnumeric);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('cutoffs',250, @isnumeric);
ip.addParameter('distType','angle');
ip.addParameter('plotHandle',[]);
ip.parse(varargin{:});
p=ip.Results;
cutoff=p.cutoffs;
kinRange=p.kinRange;

printAll=p.printAll;

%% Loading kin and EB3
%poleRefProcessIdx=(cellfun(@(p) strcmp(p.name_,'addSpindleRef'),MD.processes_));
if(~isempty(p.processSpindleRef)&&isempty(p.kinTracksWithSpindle))
    disp('reading process');tic;
    poleRefProcess=p.processSpindleRef;% MD.getProcess(find(poleRefProcessIdx,1,'last'));
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

if(isempty(kinRange))
    kinRange=1:length(kinTracksInliers);
end

%% Functionalize and processize
if(~isempty(p.processRandKin)&&isempty(p.randKinTracksWithSpindle))
    tmp=load(p.processRandKin.outFilePaths_{1}); 
    randKinTracks=tmp.randKinTracks;
else
    randKinTracks=p.randKinTracksWithSpindle;
end
%%

inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracks));

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
if(p.printAll)
    dummyProc=ExternalProcess(MD);
    printAllMTTipsKinPoleRef(kinTracksInliers(kinRange),'name','trackBased','process',dummyProc)
end;

process=p.process;
if(~isempty(process))
    procFolder=[process.getOwner().outputDirectory_  filesep 'Kin' filesep 'MTTipsBias' filesep];
    mkdir(procFolder);
    kinTracks=kinTracksInliers;
    save([procFolder 'MTTipsBiasForKin.mat'],'kinTracks');
    process.setOutFilePaths({[procFolder 'MTTipsBiasForKin.mat'],[procFolder 'MTTipsBiasForKinAndRand.mat']});
end;

if(~isempty(randKinTracks))
    randKinTracksInliers=randKinTracks(inliersKin);
    mapMTTipsToKin(randKinTracksInliers(kinRange),EB3TracksInliers,cutoff)
    appearingMTRefProjection(randKinTracksInliers(kinRange),'associatedTipsP1','associatedTipsP2');
    %    save([procFolder 'MTTipsBiasForKinAndRand.mat'],'kinTracksInliers','randKinTracksInliers');
    
    if(p.printAll)
        dummyProc=ExternalProcess(MD);
        printAllMTTipsKinPoleRef(randKinTracksInliers(kinRange),'name','randTrackBased','process',dummyProc)
    end
end

% Display for test 
% %kIdx=1;ktInlier=kinTracksInliers(kIdx);rktInlier=randKinTracksInliers(kIdx);
% 
% %[Hs]=setupFigure(1,2,2,'AxesWidth',8,'AxesHeight',4);
% %printMTTipsPoleRef(ktInlier.KP1,ktInlier.associatedTipsP2KinRef,ktInlier.associatedTipsP2Idx,'handle',Hs(1))
% %printMTTipsPoleRef(rktInlier.KP1,rktInlier.associatedTipsP2KinRef,rktInlier.associatedTipsP2Idx,'handle',Hs(2))
% 
% %%
% displayTipsBiasStat({kinTracksInliers(kinRange),randKinTracksInliers(kinRange)},{'kin','rand'});

