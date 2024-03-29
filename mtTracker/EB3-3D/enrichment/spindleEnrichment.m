function spindleEnrichment(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.addParameter('trackKTPackage',MD.searchPackageName('trackKT','selectIdx','last'));
ip.addParameter('startTime',0);
ip.addParameter('dynROIProject',[]);
ip.addParameter('useKT',true);
ip.addParameter('polarAngle',pi/6);
ip.addParameter('mappingDistance',5);
ip.addParameter('KTMoviePerc',0.5)
ip.addParameter('name',[]);
ip.addParameter('packPID',400);
ip.parse(varargin{:});
p=ip.Results;

% Process type placeholdes
packPID=p.packPID;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

packPIDTMP=packPID+1;

MD.setPackage(packPIDTMP,GenericPackage({ ...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'), ...
    ExternalProcess(MD,'project1D'),...
    ExternalProcess(MD,'project1D'),...
    ExternalProcess(MD,'spindleEnrichment'),...
    ProjectDynROIProcess(MD),...
    ProjectDynROIProcess(MD)...
    ProjectDynROIProcess(MD),...
    ProjectDynROIProcess(MD),...
    ExternalProcess(MD,'KTDetection'),...
    ExternalProcess(MD,'DetectionPolar'),...
    }));

lpid=1;

if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
   processDetectEB3=p.package.getProcess(lpid);
else
    disp('Computing  detections')
    tic
    processDetectEB3=PointSourceDetectionProcess3D(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultDetectionParams(MD,[MD.outputDirectory_ filesep 'EB3']));
    MD.addProcess(processDetectEB3);
    funParams = processDetectEB3.funParams_;
    funParams.showAll=true;
    %funParams.frameRange=[1 2];    
    funParams.alpha=0.05;
    funParams.filterSigma=[1.4396 1.4396;1.2913 1.2913  ];
    funParams.WindowSize={[],[]};
    funParams.algorithmType= {'pointSourceLM'  'pointSourceLM'};
    funParams.ConfRadius={[],[]};       
    processDetectEB3.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.isoCoord=true;
    processDetectEB3.run(paramsIn);
    toc;
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectEB3);

lpid=2;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processDetectPoles=p.package.getProcess(lpid);
else
    disp('Computing spindle ref');tic;
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
    toc;
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectPoles);

lpid=3;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processBuildRef=p.package.getProcess(lpid);
else
    processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));
    processBuildRef.run();
end
MD.getPackage(packPIDTMP).setProcess(lpid,processBuildRef);

disp('Loading spindle ref');tic;

tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);

refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;
toc;

disp('Loading and mapping +tips reference');tic;
tmp=load(processDetectEB3.outFilePaths_{1}); detection=tmp.movieInfo;

%EB3Inliers=mapDetectionsTo1DManifold(ROIs{1,2},detection,0,'distType','vertexDistOtsu');
EB3Inliers=mapDetectionsTo1DManifold(ROIs{1,2},detection,80,'distType','euclideanDist');
toc

%amiraWriteMovieInfo([fileparts(processBuildRef.outFilePaths_{2}) filesep 'AmiraDetect' filesep 'detectionLabRef.am'],EB3Inliers);

disp('Setting frame of reference');tic;
allDetections=Detections(detection);
oDetections=EB3Inliers.copy();
oDetectionsP1P2=refs(1,2).applyBase(oDetections,'');
oDetectionsP2P1=refs(2,1).applyBase(oDetections,'');
oDetectionsP1P2.addSphericalCoord();
oDetectionsP2P1.addSphericalCoord();
elevations=arrayfun(@(d,D) min(d.elevation,D.elevation),oDetectionsP1P2,oDetectionsP2P1,'unif',0);
azimuths=arrayfun(@(d,D) d.azimuth,oDetectionsP1P2,'unif',0);
toc;

lpid=4;
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processProj=p.package.getProcess(lpid);
else
processProj=ExternalProcess(MD,'project1D');
project1D(  MD, ...
            'name','fullMIPLabFrame','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',20,'fringeWidth',60, ...
            'processSingleProj',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);              
end    
MD.getPackage(packPIDTMP).setProcess(lpid,processProj);


lpid=lpid+1;
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processProjSpindleRef=p.package.getProcess(lpid);
else
    processProjSpindleRef=ExternalProcess(MD,'project1D');
    project1D(  MD,[P1,P2],'FoF',refs(1,2),'dynPoligonREF',refs(1,2).applyBase([P1,P2],''), ...
        'name','CroppedSpindleRef','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',90, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);

%% Elevation and densities

%% Debug inlier selection
% process=ExternalProcess(MD,'Detection');
% overlayProjDetectionMovie(processProj,'detections',allDetections, ...
%         'process',process,'cumulative',true,'processFrames',1, ...
%         'name','cumulDetectionDensity');
% img1=imread(sprintfPath(process.outFilePaths_{2},1));
% process=ExternalProcess(MD,'Detection');
% overlayProjDetectionMovie(processProj,'detections',EB3Inliers, ...
%         'process',process,'cumulative',true,'processFrames',1, ...
%         'name','cumulDetectionDensity');
% img2=imread(sprintfPath(process.outFilePaths_{2},1));
% figure(); imshow([img1 img2]);

%% Compute tracks if necessary
trackKTPack=p.trackKTPackage;
if(isempty(trackKTPack))
    trackKT(MD);
    trackKTPack=MD.searchPackageName('trackKT','selectIdx','last');
end

%% load and map tracks to the spindle
KTTracks=TracksHandle(trackKTPack.getProcess(3).loadChannelOutput(2));
KTTracksInliers=mapTracksTo1DManifold(ROIs{1,2},KTTracks,90,'position','start','distType','euclideanDist');

%% Filter arbitraly to weed out some weird outlies
disp(['lifetime cut off: ' num2str(MD.nFrames_*p.KTMoviePerc)]);
FilteredKTTracks=KTTracksInliers([KTTracksInliers.lifetime]>MD.nFrames_*p.KTMoviePerc);
fillTrackGaps(FilteredKTTracks);

%% Remap in Frame of reference of interest, measure spindle coordinate 
KTDetection=Detections(FilteredKTTracks.getMovieInfo());
KTDetectionsP1P2=refs(1,2).applyBase(KTDetection,'');
KTDetectionsP2P1=refs(2,1).applyBase(KTDetection,'');
KTDetectionsP1P2.addSphericalCoord();
KTDetectionsP2P1.addSphericalCoord();

%% Compute closest pole at KT end and get the associated detections for display
KTEndsIdx=arrayfun(@(t) t.tracksFeatIndxCG(end),FilteredKTTracks);
KTEndsFrame=arrayfun(@(t) t.f(end),FilteredKTTracks);
KTEndsClosestPoles=arrayfun(@(idx,f)  (KTDetectionsP1P2(f).rho(idx)>KTDetectionsP2P1(f).rho(idx))+1,KTEndsIdx,KTEndsFrame);

KTDetectionsEndIdx=arrayfun(@(f) KTEndsIdx(KTEndsFrame==f),1:length(KTDetection),'unif',0);
KTEnds=KTDetection.getSelectIdx(KTDetectionsEndIdx);

%% Strips unlinked detection and recompute projection (for integrity)
linkedDet=arrayfun(@(d) ~isnan(d.xCoord(:,1)),KTDetection,'unif',0);
KTDetection.selectIdx(linkedDet);
KTDetectionsP1P2=refs(1,2).applyBase(KTDetection,'');
KTDetectionsP2P1=refs(2,1).applyBase(KTDetection,'');
KTDetectionsP1P2.addSphericalCoord();
KTDetectionsP2P1.addSphericalCoord();

KTelevations=arrayfun(@(d,D) min(d.elevation,D.elevation),KTDetectionsP1P2,KTDetectionsP2P1,'unif',0);
KTpoleDistances=arrayfun(@(d,D) min(d.rho,D.rho),KTDetectionsP1P2,KTDetectionsP2P1,'unif',0);
KTazimuths=arrayfun(@(d,D) d.azimuth,KTDetectionsP1P2,'unif',0);
KTTimesSteps=arrayfun(@(d,t) p.startTime+MD.timeInterval_*(t-1)*ones(size(d.xCoord(:,1))),KTDetection,1:length(KTDetection),'unif',0);

%%
[P1KTROI]=build1DManifold([P1],FilteredKTTracks(KTEndsClosestPoles==1));
[P2KTROI]=build1DManifold([P2],FilteredKTTracks(KTEndsClosestPoles==2));
PKTROI=[P1KTROI P2KTROI];
disp('Mapping +Tips on PTK manifolds');tic;
PKTMapDist=p.mappingDistance;
% Map each dyn ROI and measure distance
allKTMappedTips=cell(1,numel(PKTROI));
allKTMappedTipsDist=cell(1,numel(PKTROI));
parfor mIdx=1:numel(PKTROI);
    [~,indices,distance]=mapDetectionsTo1DManifold(PKTROI{mIdx},oDetections,PKTMapDist, ... 
                                                        'distType','normalDistPseudOptimized');
    allKTMappedTips{mIdx}=indices;
    allKTMappedTipsDist{mIdx}=distance;
end

%% Compute the Union of all the mapped +Tips (sub-optimal, should be done per detection with a quadtree)
KTMappedTips=allKTMappedTips{1};
KTMappedIdx=cellfun(@(i) double(i),allKTMappedTips{1},'unif',0);
KTMinDist=allKTMappedTipsDist{1};
for mIdx=2:numel(PKTROI);
%    KTMappedTips=cellfun(@(m,i) (uint16(m)+uint16(i)),KTMappedTips,allKTMappedTips{mIdx},'unif',0);
    KTMappedTips=cellfun(@(m,i) (m|i),KTMappedTips,allKTMappedTips{mIdx},'unif',0);
    KTMinDist=cellfun(@(m,d) min(m,d),KTMinDist,allKTMappedTipsDist{mIdx},'unif',0);
   %KTMappedIdx=cellfun(@(i,m,d) i.*(1-(m==d))+(m==d)&(m>0)*mIdx,KTMappedIdx,KTMinDist,allKTMappedTipsDist{mIdx},'unif',0);
    KTMappedIdx=cellfun(@(i,m,d)  i+double(m&(i==0))*mIdx,KTMappedIdx,KTMappedTips,allKTMappedTipsDist{mIdx},'unif',0);

end
figure;histogram(KTMappedTips{2})
figure;histogram(KTMinDist{2})
figure;histogram(allKTMappedTipsDist{2}{2})
polarIndices=KTMappedTips;
toc;

% For each dynROI estimate density

%% Create inverse PKT manifolds and map detection to them
% Mirroring KTs
invFilteredKTTracks=FilteredKTTracks.getAddCoord(P1.getMultCoord(-1)).multCoord(-1).addCoord(P1);
[InvP1KTROI]=build1DManifold([P1],invFilteredKTTracks);
invFilteredKTTracks=FilteredKTTracks.getAddCoord(P2.getMultCoord(-1)).multCoord(-1).addCoord(P2);
[InvP2KTROI]=build1DManifold([P2],invFilteredKTTracks);
PKTROIAstral=[InvP1KTROI InvP2KTROI];
disp('Mapping +Tips on inversed PTK manifolds');
tic;
PKTMapDist=p.mappingDistance;
allKTMappedTips=cell(1,numel(PKTROIAstral));
parfor mIdx=1:numel(PKTROIAstral);
    [~,indices]=mapDetectionsTo1DManifold(PKTROIAstral{mIdx},oDetections,PKTMapDist, ... 
                                                        'distType','normalDistPseudOptimized');
    allKTMappedTips{mIdx}=indices;
end
% Compute the Union of all the mapped +Tips
KTMappedTips=allKTMappedTips{1};
for mIdx=2:numel(PKTROIAstral);
    KTMappedTips=cellfun(@(m,i) (m|i),KTMappedTips,allKTMappedTips{mIdx},'unif',0);
end
astralIndices=KTMappedTips;
toc;

process_PPKSlice=MD.findProcessTag('project1D_PPKSlice','safeCall',true);
if(isempty(process_PPKSlice))
    process_PPKSlice=ProjectDynROIProcess(MD,'PPKTSlice');
    KT=FilteredKTTracks(find([FilteredKTTracks.lifetime]==MD.nFrames_,1));
    refPPKT=copy(refs(1,2));
    refPPKT.genBaseFromZ(KT);
    KTRefPPKT=refPPKT.applyBase(KT,'');
    KTRefPPKTOpposite=KTRefPPKT.copy();
    KTRefPPKTOpposite.x=-KTRefPPKTOpposite.x;
    project1D(  MD,[],'FoF',refPPKT,'dynPoligonREF',[refPPKT.applyBase([P1,P2],'') KTRefPPKT   KTRefPPKTOpposite], ...
        'name','PPKSlice','channelRender','grayRed','saveSingleProj',true, 'fringeWidth',[40,5,40], ...
        'processSingleProj',process_PPKSlice, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
    process_PPKSlice.setProcessTag(['project1D_PPKSlice']);
    MD.addProcess(process_PPKSlice);
end

% process_PPKSlice=MD.findProcessTag('projDynROI_PPKSlice','safeCall',true);
% if(isempty(process_PPKSlice))
%     % Building ROI
%     KT=FilteredKTTracks(find([FilteredKTTracks.lifetime]==MD.nFrames_,1));
%     refPPKT=copy(refs(1,2));
%     refPPKT.genBaseFromZ(KT);
%     KTRefPPKT=refPPKT.applyBase(KT,'');
%     KTRefPPKTOpposite=KTRefPPKT.copy();
%     KTRefPPKTOpposite.x=-KTRefPPKTOpposite.x;

%     % 


%     process_PPKSlice = ProjectDynROIProcess(MD,'projDynROI_PPKSlice');


%     projectDynROI(  MD,[],[refPPKT.applyBase([P1,P2],'') KTRefPPKT   KTRefPPKTOpposite], ...
%         'FoF',refPPKT,'name','PPKSlice','channelRender','grayRed','saveSingleProj',true, 'fringeWidth',[40,5,40], ...
%         'processSingleProj',process_PPKSlice, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
%     process_PPKSlice.setProcessTag(['projDynROI_PPKSlice']);
%     MD.addProcess(process_PPKSlice);
% end

% detection debug visualization
%process=ExternalProcess(MD,'Detection polar');
% overlayProjDetectionMovie(processProjSpindleRef,'detections',oDetectionsP1P2.getSelectIdx(KTMappedTips),'process',process, ... 
%             'name','KTMappedTips');
% img=imread(sprintfPath(process.outFilePaths_{4},1));
% figure();imshow([img]);

disp('Compute densities in the spindle'); tic;
densities=estimateDensity(oDetections,20);
poleDistances=arrayfun(@(d,D) min(d.rho,D.rho),oDetectionsP1P2,oDetectionsP2P1,'unif',0);

%% Define astral and polar region in a naive way
% astralThresh=-pi/4;
% polarThresh=pi/4;
astralThresh=-p.polarAngle;
polarThresh=p.polarAngle;

%% Count and density function of distance to poles
% astralDistBinning=linspace(5,50,40);
% polarDistBinning=linspace(5,50,40);
astralDistBinning=5:1:50;
polarDistBinning=5:1:50;

astralDistBinning=1:1:100;
polarDistBinning=1:1:100;

inDistBound=cellfun(@(p)  ((p<max(polarDistBinning))&(p>min(polarDistBinning)))' ,poleDistances,'unif',0);

polarIndices=cellfun(@(p,e)  p&e ,inDistBound,polarIndices,'unif',0);
astralIndices=cellfun(@(p,e)  p&e ,inDistBound,astralIndices,'unif',0);

polarIndicesAngle=cellfun(@(p,e)  p&(e'>polarThresh) ,inDistBound,elevations,'unif',0);
astralIndicesAngle=cellfun(@(p,e)  p&(e'<astralThresh) ,inDistBound,elevations,'unif',0);

OutsideKTPolar=cellfun(@(k,a) (~k&a),polarIndices,polarIndicesAngle,'unif',0);
OutsideKTAstral=cellfun(@(k,a) (~k&a),astralIndices,astralIndicesAngle,'unif',0);

if(~p.useKT)
  polarIndices=polarIndicesAngle;
  astralIndices=astralIndicesAngle;
end
toc;

disp('Order  stats in the spindle');tic;
allElevs=vertcat(elevations{:});
allDens=vertcat(densities{:});
amps=arrayfun(@(d) d.amp,oDetections,'unif',0); % just for data structure coherence

scoresBin=-pi/2:0.05:pi/2;

%% Density vs Elevation stats
[counts,edges,binIdx]=histcounts(allElevs,scoresBin);
[means,meds,stds,orderedIndex,counts] = statPerIndx(allDens,binIdx+1,1:length(scoresBin));

%% total count and pole distance/elevation

allDist=vertcat(poleDistances{:});
allAmp=vertcat(amps{:});
allAmp=allAmp(:,1);
timesSteps=arrayfun(@(d,t) p.startTime+MD.timeInterval_*(t-1)*ones(size(d.xCoord(:,1))),oDetections,1:length(oDetections),'unif',0);

polarDetectionsMask=logical(horzcat(polarIndices{:}));
astralDetectionsMask=logical(horzcat(astralIndices{:}));

[polarCounts,polarEdges,polarBinIdx]=histcounts(allDist(polarDetectionsMask),polarDistBinning);
[meanPolarDensities,medPolarDensities,stdPolarDensities]=statPerIndx(allDens(polarDetectionsMask),polarBinIdx,1:length(polarDistBinning)-1);
[meanPolarAmp,medPolarAmp,stdPolarAmp,~,polarCount,sumPolarAmp]=statPerIndx(allAmp(polarDetectionsMask),polarBinIdx,1:length(polarDistBinning)-1);
coneVolumes=2*pi*(1-cos(pi/2-polarThresh))*polarDistBinning.^3/3;
polarVolumes=(coneVolumes(2:end)-coneVolumes(1:(end-1)));
normalizedPolarAmp=sumPolarAmp./polarVolumes;
normalizedPolarCount=polarCounts./polarVolumes;

[astralCounts,astralEdges,astralBinIdx]=histcounts(allDist(astralDetectionsMask),astralDistBinning);
[meanAstralDensities,medAstralDensities,stdAstralDensities]=statPerIndx(allDens(astralDetectionsMask),astralBinIdx,1:length(astralDistBinning)-1);
[meanAstralAmp,medAstralAmp,stdAstralAmp,~,astralCount,sumAstralAmp]=statPerIndx(allAmp(astralDetectionsMask),astralBinIdx,1:length(astralDistBinning)-1);
coneVolumes=2*pi*(1-cos(pi/2+astralThresh))*astralDistBinning.^3/3;
astralVolumes=(coneVolumes(2:end)-coneVolumes(1:(end-1)));
normalizedAstralAmp=sumAstralAmp./astralVolumes;
normalizedAstralCount=astralCounts./astralVolumes;


%% Density function of distance to poles
densityLimit=[0 max([meanAstralDensities(:); meanPolarDensities(:)])*1.05];
toc;

disp('Save final process');tic;
lpid=lpid+1;
processEnrich=ExternalProcess(MD,'spindleEnrichment');
outputFileName=[MD.outputDirectory_ filesep 'enrichment' filesep p.name filesep 'elevation-density-poleDistance.mat'];
mkdirRobust(fileparts(outputFileName));
save(outputFileName,'elevations','azimuths','densities','timesSteps','poleDistances','amps', ... 
    'astralThresh','polarThresh','polarDistBinning','astralDistBinning','densityLimit', ...
    'KTelevations','KTpoleDistances','KTazimuths','KTTimesSteps','astralDetectionsMask','polarDetectionsMask');
outputDirPlot=[MD.outputDirectory_ filesep 'enrichment' filesep p.name filesep 'plot' filesep];
processEnrich.setOutFilePaths({ ... 
    [outputDirPlot  'ElevationDens.png'], ...
     outputFileName, ...
    [outputDirPlot  'poleDist.png'], ...
    [outputDirPlot  'densityVsDist.png'], ...
    [outputDirPlot  'ampVsDist.png'], ...
    [outputDirPlot  'normAmpVsDist.png'], ...
    [outputDirPlot  'normCountVsDist.png'], ...
    });
toc;
MD.getPackage(packPIDTMP).setProcess(lpid,processEnrich);

if(GenericPackage.processExist(p.package,lpid))
    processEnrich=p.package.getProcess(lpid);
    %load(processEnrich.outFilePaths_{2});
else
    disp('Refresh stats display');tic;
    % [Handle,~,F]=setupFigure(1,1,1);
    % shadedErrorBar(scoresBin,means,stds,'r',1);       
    % xlabel('Elevation (rad)');
    % xlim(minmax(scoresBin));
    % ylim([0,8]);
    % ylabel('Comet density (1/\mu m^3)')
    % mkdirRobust(outputDirPlot);
    % print([outputDirPlot  'ElevationDens.png'],'-dpng');
    % print([outputDirPlot  'ElevationDens.eps'],'-depsc');

    % Count function of distance to poles
    [Handle,~,F]=setupFigure(1,2,2);
    plot(Handle(1),astralEdges(1:end-1)*0.1,astralCounts/MD.nFrames_);
    plot(Handle(2),polarEdges(1:end-1)*0.1,polarCounts/MD.nFrames_);
    arrayfun(@(h) ylim(h,[0 12]),Handle);
    % arrayfun(@(h) ylim(h,[0 max([astralCounts(:); polarCounts(:)]/MD.nFrames_)*1.05]),Handle);
    arrayfun(@(h) xlim(h,[1,6.5]),Handle);
    xlabel(Handle(1),'Pole distance (\mum)');
    xlabel(Handle(2),'Pole distance (\mum)');
    ylabel(Handle(1),'Astral Comet count')
    ylabel(Handle(2),'Polar Comet count')
    mkdirRobust(outputDirPlot);
    print([outputDirPlot  'poleDist.png'],'-dpng');
    print([outputDirPlot  'poleDist.eps'],'-depsc');


%     [Handle,~,F]=setupFigure(1,2,2);
%     plot(Handle(1),astralEdges(1:end-1)*0.1,meanAstralDensities);
%     plot(Handle(2),polarEdges(1:end-1)*0.1,meanPolarDensities);
%     arrayfun(@(h) xlim(h,[0,6.5]),Handle);
%     xlabel(Handle(1),'Pole distance (\mum)');
%     xlabel(Handle(2),'Pole distance (\mum)');
%     ylabel(Handle(1),'Astral Comet density (per volume)')
%     ylabel(Handle(2),'Polar Comet density (per volume)')
%     arrayfun(@(h) ylim(h,densityLimit),Handle);
%     mkdirRobust(outputDirPlot);
%     printPNGEPSFIG(F,outputDirPlot,'densityVsDist');

% % Amplitude function of distance to poles
% [Handle,~,F]=setupFigure(1,2,2);
% plot(Handle(1),astralEdges(1:end-1)*0.1,meanAstralAmp);
% plot(Handle(2),polarEdges(1:end-1)*0.1,meanPolarAmp);
% arrayfun(@(h) xlim(h,[0,6.5]),Handle);
% xlabel(Handle(1),'Pole distance (\mum)');
% xlabel(Handle(2),'Pole distance (\mum)');
% ylabel(Handle(1),'Astral Comet Mean Intensity')
% ylabel(Handle(2),'Polar Comet Mean Intensity')
% arrayfun(@(h) ylim(h,[min([meanAstralAmp(:); meanPolarAmp(:)])*0.95 max([meanAstralAmp(:); meanPolarAmp(:)])*1.05]),Handle);
% printPNGEPSFIG(F,outputDirPlot,'ampVsDist');

% %% Normalized amplitude function of distance to poles
% [Handle,~,F]=setupFigure(1,2,2);
% plot(Handle(1),astralEdges(1:end-1)*0.1,normalizedAstralAmp);
% plot(Handle(2),polarEdges(1:end-1)*0.1,normalizedPolarAmp);
% arrayfun(@(h) xlim(h,[1,6.5]),Handle);
% xlabel(Handle(1),'Pole distance (\mum)');
% xlabel(Handle(2),'Pole distance (\mum)');
% ylabel(Handle(1),'Astral Comet Intensity (Solid Angle Norm)')
% ylabel(Handle(2),'Polar Comet Intensity (Solid Angle Norm)')
% arrayfun(@(h) ylim(h,[min([normalizedPolarAmp(:); normalizedAstralAmp(:)])*0.95 max([normalizedPolarAmp(:); normalizedAstralAmp(:)])*1.05]),Handle);
% printPNGEPSFIG(F,outputDirPlot,'normAmpVsDist');

% %% Normalized count function of distance to poles
% [Handle,~,F]=setupFigure(1,2,2);
% plot(Handle(1),astralEdges(1:end-1)*0.1,normalizedAstralCount);
% plot(Handle(2),polarEdges(1:end-1)*0.1,normalizedPolarCount);
% arrayfun(@(h) xlim(h,[1,6.5]),Handle);
% xlabel(Handle(1),'Pole distance (\mum)');
% xlabel(Handle(2),'Pole distance (\mum)');
% ylabel(Handle(1),'Astral Comet density (normCount)')
% ylabel(Handle(2),'Polar Comet density (normCount)')
% arrayfun(@(h) ylim(h,[0 max([normalizedPolarCount(:); normalizedAstralCount(:)])*1.05]),Handle);
% printPNGEPSFIG(F,outputDirPlot,'normCountVsDist');

% %% volumes function of distance to poles
% [Handle,~,F]=setupFigure(1,2,2);
% plot(Handle(1),astralEdges(1:end-1)*0.1,astralVolumes);
% plot(Handle(2),polarEdges(1:end-1)*0.1,polarVolumes);
% arrayfun(@(h) xlim(h,[1,6.5]),Handle);
% xlabel(Handle(1),'Pole distance (\mum)');
% xlabel(Handle(2),'Pole distance (\mum)');
% ylabel(Handle(1),'Astral area')
% ylabel(Handle(2),'Polar area')
% arrayfun(@(h) ylim(h,[0 max([astralVolumes(:); polarVolumes(:)])*1.05]),Handle);
% printPNGEPSFIG(F,outputDirPlot,'polarVolumes');
toc;
end

if(~isempty(p.dynROIProject))
    processProjDynROI=p.dynROIProject;
else
    processProjDynROI=process_PPKSlice;
end

%% Rendering cumulative distribution for debugging and illustrative purposes
myColormap=255*jet(256);
myColormap(1,:)=[200 200 200];
lpid=lpid+1;
if(~GenericPackage.processExist(p.package,lpid))
    disp('Render Astral  still picture');tic;
    process=ProjRendering(processProjDynROI,['astralCumulPoleDensity' p.name]);

    astralDetections=oDetections.copy().selectIdx(astralIndices);
    astralDensitiesIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),densityLimit))+1 , densities,astralIndices,'unif',0);
    astralDistanceIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),[20 60]))+1 , poleDistances,astralIndices,'unif',0);

    overlayProjDetectionMovie(processProjDynROI,'detections',astralDetections, 'process',process, ... 
        'cumulative',true,'processFrames',1, ... 
        'colorIndx',astralDistanceIndx,'colormap',myColormap,'name',['astralCumulPoleDensity' p.name]);
    if(p.useKT)
        processEB=copy(process);
        tracks=fillTrackGaps(FilteredKTTracks);
        KTFilteredDetection=Detections(tracks.getMovieInfo());
        overlayProjDetectionMovie(processEB,'detections',KTEnds, 'process',process, ... 
            'cumulative',true,'processFrames',1,'name','KTDetection','radius',p.mappingDistance);
    end
    figure();imshow(process.loadFrame(1,1));
    toc
else
    process=p.package.getProcess(lpid);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

lpid=lpid+1;
if(~GenericPackage.processExist(p.package,lpid))
    disp('Render Polar still picture');tic;
    process=ProjRendering(processProjDynROI,['polarCumulPoleDensity' p.name]);
    polarDetections=oDetections.copy().selectIdx(polarIndices);
    polarDensitiesIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),densityLimit))+1 , densities,polarIndices,'unif',0);
    polarDistanceIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),[20 60]))+1 , poleDistances,polarIndices,'unif',0);
    polarKTIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),[1 numel(PKTROI)]))+1 , KTMappedIdx,polarIndices,'unif',0);
    % polarKTIndx=cellfun(@(m,i) floor(254*i.*mat2gray(m,[1 numel(PKTROI)]))+1 , KTMappedIdx,polarIndices,'unif',0);

    overlayProjDetectionMovie(processProjDynROI,'detections',polarDetections, 'process',process, ... 
        'cumulative',true,'processFrames',1, ...
        'colorIndx',polarDistanceIndx,'colormap',myColormap,'name',['polarCumulPoleDensity' p.name]);
    if(p.useKT) 
        processEB=copy(process);
        overlayProjDetectionMovie(processEB,'detections',KTEnds, 'process',process, ... 
            'cumulative',true,'processFrames',1,'name','KTDetection','radius',p.mappingDistance);
    end
    figure();imshow(process.loadFrame(1,1));
    toc;
else
    process=p.package.getProcess(lpid);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);



lpid=lpid+1;
KTColormap=255*winter(256);
if(~GenericPackage.processExist(p.package,lpid))
    disp('Render Astral  detections in dynROI');tic;
    process=ProjRendering(processProjDynROI,['astralPoleDensity-' p.name]);
    % astralDetections=oDetections.copy().selectIdx(astralIndices);
    % astralDistancesIndx=cellfun(@(e,i)  floor(254*mat2gray(e(i)))+1, poleDistances,polarIndx,'unif',0);
    astralDistanceIndx=cellfun(@(e,i) floor(254*i'.*mat2gray(e,[0 60]))+1 , poleDistances,astralIndices,'unif',0);
    astralDensitiesIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),densityLimit))+1 , densities,astralIndices,'unif',0);
    overlayProjDetectionMovie(processProjDynROI,'detections',oDetections, ...
         'process',process, 'cumulative',false, ... 
        'colorIndx',astralDistanceIndx,'colormap',myColormap,'name',['astralPoleDensity-' p.name]);
    if(p.useKT)
        tracks=fillTrackGaps(FilteredKTTracks);
        processEB=copy(process);
        KTFilteredDetection=Detections(tracks.getMovieInfo());
        overlayProjDetectionMovie(processEB,'detections',KTFilteredDetection, 'process',process, ... 
            'cumulative',false,'name','KTDetection','radius',p.mappingDistance);
    end
    toc;
else
    process=p.package.getProcess(lpid);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

%%
poleDistLimit=[19.999 20];
lpid=lpid+1;
if(~GenericPackage.processExist(p.package,lpid))
    disp('Render polar  detections in dynROI');tic;
    process=ProjRendering(processProjDynROI,['polarPoleDensity-' p.name]);
    polarDetections=oDetections.copy().selectIdx(polarIndices);
    % poleDistLimit used to validate detection.

    polarDistanceIndx=cellfun(@(e,i) floor(254*i'.*mat2gray(e,[0 60]))+1 , poleDistances,polarIndices,'unif',0);
    % polarDistanceIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),poleDistLimit))+1 , poleDistances,polarIndices,'unif',0);
    % polarDistanceIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),[0 60]))+1 , poleDistances,polarIndices,'unif',0);
    polarDensitiesIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),densityLimit))+1 , densities,polarIndices,'unif',0);
    polarKTIndx=cellfun(@(e,i) floor(254*i'.*mat2gray(e',[1 numel(PKTROI)]))+1 , KTMappedIdx,polarIndices,'unif',0);
    overlayProjDetectionMovie(processProjDynROI,'detections',oDetections, 'process',process, ... 
        'cumulative',false, ...
        'colorIndx',polarDistanceIndx,'colormap',myColormap,'name',['polarPoleDensity-' p.name],'radius',2);
    if(p.useKT)
        % overlayProjTracksMovie(process,'tracks', ref.applyBase(FilteredKTTracks,''), ... 
        %     'colorIndx',ceil(255*mat2gray([FilteredKTTracks.lifetime]',[1 150]))+1, ... 
        %     'dragonTail',10,'colormap',KTColormap,'name',['KT' p.name],'process',process);
        tracks=fillTrackGaps(FilteredKTTracks);
        processEB=copy(process);
        KTFilteredDetection=Detections(tracks.getMovieInfo());
        overlayProjDetectionMovie(processEB,'detections',KTFilteredDetection, 'process',process, ... 
            'cumulative',false,'name','KTDetection','radius',p.mappingDistance);
    end

    figure();imshow(process.loadFrame(1,1));
    toc;
else
    process=p.package.getProcess(lpid);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

lpid=lpid+1;
if(~GenericPackage.processExist(p.package,lpid))
    disp('Render KT detections (TrackKT is broken ...)');tic;
    process=ProjectDynROIProcess(MD);
    overlayProjDetectionMovie(processProjSpindleRef,'detections',KTDetections, 'process',process, ... 
        'cumulative',false,'name',['KTDetection-' p.name]);
    toc;
else
    process=p.package.getProcess(lpid);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

lpid=lpid+1;
if(~GenericPackage.processExist(p.package,lpid))
    disp('Rendering +Tips non-mapped by KT in the interpolar area');tic;
    process=ExternalProcess(MD,'Detection polar');
    polarDetections=oDetections.copy().selectIdx(OutsideKTPolar);
    polarDensitiesIndx=cellfun(@(e,i) floor(254*mat2gray(e(i),densityLimit))+1 , densities,OutsideKTPolar,'unif',0);
    overlayProjDetectionMovie(processProjDynROI,'detections',refs(1,2).applyBase(polarDetections,''), 'process',process, ... 
        'cumulative',true,'processFrames',1, ...
        'colorIndx',polarDensitiesIndx,'colormap',myColormap,'name',['polarCumulNonMapped' p.name]);
    if(p.useKT)

        overlayProjDetectionMovie(process,'detections',refs(1,2).applyBase(KTEnds,''), 'process',process, ... 
            'cumulative',true,'processFrames',1 ...
            ,'name','CumulPoleDensityKT');
    end
    imgPolar=imread(sprintfPath(process.outFilePaths_{4},1));
    figure();
    imshow([imgPolar]);
    toc;
else
    process=p.package.getProcess(lpid);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

MD.setPackage(packPID,MD.getPackage(packPIDTMP));

function [means,meds,stds,orderedIndex,counts,sums] = statPerIndx(values,index,indexValues)
means=zeros(1,max(index));
meds=zeros(1,max(index));
stds=zeros(1,max(index));
counts=zeros(1,max(index));
sums=zeros(1,max(index));

orderedIndex=indexValues;
for t=indexValues
    if(t>0)
    biasAtTime=(values(index==t));
    means(t)=mean(biasAtTime);
    meds(t)=median(biasAtTime);
    stds(t)=std(biasAtTime);
    counts(t)=sum(index==t);
    sums(t)=sum(biasAtTime);
    end
end