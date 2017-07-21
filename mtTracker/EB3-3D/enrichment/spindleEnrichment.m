function spindleEnrichment(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.parse(varargin{:});
p=ip.Results;

% Process type placeholdes
packPID=400;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

packPIDTMP=packPID+1; MD.setPackage(packPIDTMP,GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'), ...
    ExternalProcess(MD,'project1D'),...
    ExternalProcess(MD,'project1D'),...    
    ExternalProcess(MD,'spindleEnrichment')
  }));

lpid=1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
   processDetectEB3=p.package.getProcess(lpid);
else
    processDetectEB3=PointSourceDetectionProcess3D(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultDetectionParams(MD,[MD.outputDirectory_ filesep 'EB3']));
    MD.addProcess(processDetectEB3);
    funParams = processDetectEB3.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.filterSigma=[1.4396 1.4396;1.2913 1.2913  ];
    funParams.WindowSize={[],[]};
    funParams.algorithmType= {'pointSourceLM'  'pointSourceLM'};
    funParams.ConfRadius={[],[]};       
    processDetectEB3.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.isoCoord=true;
    processDetectEB3.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectEB3);

lpid=2;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processDetectPoles=p.package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
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

tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);

refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;

tmp=load(processDetectEB3.outFilePaths_{1}); detection=tmp.movieInfo;

EB3Inliers=mapDetectionsTo1DManifold(ROIs{1,2},detection,0,'distType','vertexDistOtsu');

amiraWriteMovieInfo([fileparts(processBuildRef.outFilePaths_{2}) filesep 'AmiraDetect' filesep 'detectionLabRef.am'],EB3Inliers);

oDetections=Detections(EB3Inliers);
oDetectionsP1P2=refs(1,2).applyBase(oDetections,'');
oDetectionsP2P1=refs(2,1).applyBase(oDetections,'');
oDetectionsP1P2.addSphericalCoord();
oDetectionsP2P1.addSphericalCoord();
elevations=arrayfun(@(d,D) min(d.elevation,D.elevation),oDetectionsP1P2,oDetectionsP2P1,'unif',0);

lpid=4;
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processProj=p.package.getProcess(lpid);
else
processProj=ExternalProcess(MD,'rawProj');
project1D(  MD, ...
            'name','fullMIPLabFrame','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',20,'fringeWidth',60, ...
            'processSingleProj',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);              
end    
MD.getPackage(packPIDTMP).setProcess(lpid,processProj);


lpid=lpid+1;
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processProjSpindleRef=p.package.getProcess(lpid);
else
    processProjSpindleRef=ExternalProcess(MD,'rawProj');
    project1D(  MD,[P1,P2],'FoF',refs(1,2),'dynPoligonREF',refs(1,2).applyBase([P1,P2],''), ...
        'name','CroppedSpindleRef','channelRender','grayRed','saveSingleProj',true, 'insetFringeWidth',80, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);

%% Elevation and densities
myColormap=255*jet(256);
densities=estimateDensity(oDetections,30);

% overlayProjDetectionMovie(processProj,'detections',oDetections, ... 
%             'colorIndx',cellfun(@(e) 1+ceil(mat2gray(e)*255), elevations,'unif',0),'colormap',myColormap,'name','elevations')
% 
% overlayProjDetectionMovie(processProj,'detections',oDetections, ... 
%             'colorIndx',cellfun(@(e) 1+ceil(mat2gray(e)*255), densities,'unif',0),'colormap',myColormap,'name','densities')

allElevs=vertcat(elevations{:});
allDens=vertcat(densities{:});
%
scoresBin=-pi/2:0.05:pi/2;
[counts,edges,binIdx]=histcounts(allElevs,scoresBin);
[means,meds,stds,orderedIndex,counts] = statPerIndx(allDens,binIdx+1);
% [H,~,F]=setupFigure(1,2,2);
% scatter(H(1),vertcat(elevations{:}),vertcat(densities{:}));
% axes(H(2));
% shadedErrorBar(scoresBin,means,stds,'r',1);       

[Handle,~,F]=setupFigure(1,1,1);
shadedErrorBar(scoresBin,means,stds,'r',1);       
xlabel('Elevation (rad)');
xlim(minmax(scoresBin));
ylim([0,8]);

ylabel('Comet density (1/\mu m^3)')
outputDirPlot=[MD.outputDirectory_ filesep 'enrichment' filesep 'plot' filesep];
mkdirRobust(outputDirPlot);
print([outputDirPlot  'ElevationDens.png'],'-dpng');
print([outputDirPlot  'ElevationDens.eps'],'-depsc');

%% total count and pole distance/elevation

poleDistances=arrayfun(@(d,D) min(d.rho,D.rho),oDetectionsP1P2,oDetectionsP2P1,'unif',0);
allDist=vertcat(poleDistances{:});
%%
% myColormap=255*jet(256);
% overlayProjDetectionMovie(processProj,'detections',oDetections, ... 
%             'colorIndx',cellfun(@(p,e) 1+ceil((e<pi/6).*mat2gray(p)*255), poleDistances,elevations,'unif',0),'colormap',myColormap,'name','AstralPoleDistances')
% overlayProjDetectionMovie(processProj,'detections',oDetections, ... 
%             'colorIndx',cellfun(@(p,e) 1+ceil((e>pi/6).*mat2gray(p)*255), poleDistances,elevations,'unif',0),'colormap',myColormap,'name','polarPoleDistances')        
%%
[polarCounts,polarEdges,binIdx]=histcounts(allDist(allElevs>pi/6),100);
[astralCounts,astralEdges,binIdx]=histcounts(allDist(allElevs<pi/6),100);
   

[Handle,~,F]=setupFigure(1,2,2);
plot(Handle(1),astralEdges(1:end-1)*0.1,astralCounts);
plot(Handle(2),polarEdges(1:end-1)*0.1,polarCounts);
xlabel(Handle(1),'Pole distance (\mum)');
xlabel(Handle(2),'Pole distance (\mum)');

ylabel(Handle(1),'Astral Comet count')
ylabel(Handle(2),'Polar Comet count')
mkdirRobust(outputDirPlot);
print([outputDirPlot  'poleDist.png'],'-dpng');
print([outputDirPlot  'poleDist.eps'],'-depsc');

lpid=lpid+1;
processEnrichVsElev=ExternalProcess(MD,'spindleEnrichment');
outputFileName=[MD.outputDirectory_ filesep 'enrichment' filesep 'elevation-density-poleDistance.mat'];
save(outputFileName,'elevations','densities','poleDistances','Handle')
processEnrichVsElev.setOutFilePaths({[outputDirPlot  'ElevationDens.png'],outputFileName,[outputDirPlot  'poleDist.png']});
MD.getPackage(packPIDTMP).setProcess(lpid,processEnrichVsElev);

MD.setPackage(packPID,MD.getPackage(packPIDTMP));

function [means,meds,stds,orderedIndex,counts] = statPerIndx(values,index)
means=zeros(1,max(index));
meds=zeros(1,max(index));
stds=zeros(1,max(index));
counts=zeros(1,max(index));

orderedIndex=1:max(index);
for t=unique(index)'
    biasAtTime=(values(index==t));
    means(t)=mean(biasAtTime);
    meds(t)=median(biasAtTime);
    stds(t)=std(biasAtTime);
    counts(t)=sum(index==t);
end