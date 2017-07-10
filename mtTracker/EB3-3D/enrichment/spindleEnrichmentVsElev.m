function spindleEnrichmentVsElev(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.parse(varargin{:});
p=ip.Results;

% Process type placeholdes
packPID=400;
MD.setPackage(packPID,GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'spindleEnrichmentVsElev'), ...
    ExternalProcess(MD,'project1D'),...
    ExternalProcess(MD,'spindleEnrichmentVsElev')
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
MD.getPackage(packPID).setProcess(lpid,processDetectEB3);

lpid=2;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processDetectPoles=p.package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end
MD.getPackage(packPID).setProcess(lpid,processDetectPoles);

lpid=3;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processBuildRef=p.package.getProcess(lpid);
else
    processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));
    processBuildRef.run();
end
MD.getPackage(packPID).setProcess(lpid,processBuildRef);

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
MD.getPackage(packPID).setProcess(lpid,processProj);

myColormap=255*jet(256);
overlayProjDetectionMovie(processProj,'detections',oDetections, ... 
            'colorIndx',cellfun(@(e) 1+ceil(mat2gray(e)*255), elevations,'unif',0),'colormap',myColormap)
densities=estimateDensity(oDetections,30);
overlayProjDetectionMovie(processProj,'detections',oDetections, ... 
            'colorIndx',cellfun(@(e) 1+ceil(mat2gray(e)*255), densities,'unif',0),'colormap',myColormap,'name','densities')

allElevs=vertcat(elevations{:});
allDens=vertcat(densities{:});
%%
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
outputDirPlot=[MD.outputDirectory_ filesep 'density' filesep 'plot' filesep];
mkdirRobust(outputDirPlot);
print([outputDirPlot  'ElevationDens.png'],'-dpng');
print([outputDirPlot  'ElevationDens.eps'],'-depsc');

lpid=5;
processEnrichVsElev=ExternalProcess(MD,'spindleEnrichmentVsElev');
save([MD.outputDirectory_ filesep 'density' filesep 'elevation-density.mat'],'elevations','densities','Handle')
processEnrichVsElev.setOutFilePaths({[MD.outputDirectory_ filesep 'density' filesep 'elevation-density.mat'],[outputDirPlot  'ElevationDens.png']});
MD.getPackage(packPID).setProcess(lpid,processEnrichVsElev);

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