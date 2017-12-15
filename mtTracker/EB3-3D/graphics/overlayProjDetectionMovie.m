function overlayProjDetectionMovie(processProj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('processProj');
ip.addOptional('detections',[]);
ip.addOptional('colormap',[]);
ip.addOptional('process',[]);
ip.addOptional('processFrames',[]);
ip.addOptional('cumulative',false);
ip.addOptional('colorIndx',[],@iscell);
ip.addOptional('name','detections');
ip.parse(varargin{:});
p=ip.Results;

detections=p.detections;
cumulative=p.cumulative;
%% testing imwarp to crop the image

projDataIdx=5;
ref=[];

if(isa(processProj,'ExternalProcess'))
  processProjDynROI=ProjectDynROIRendering();
  processProjDynROI.importFromDeprecatedExternalProcess(processProj);
  try % Handle Project1D/ProjDyn different outFilePaths_spec (need to be defined through a class...)
    projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  catch
    projDataIdx=4;
  end
  projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  processProjDynROI.setBoundingBox( [projData.minXBorder projData.maxXBorder], [projData.minYBorder projData.maxYBorder], [projData.minZBorder projData.maxZBorder]);
  processProj=processProjDynROI;
end


ref=get(processProj,'ref');
detections=ref.applyBase(detections,'');  
projData=processProj;


frameNb=min([projData.frameNb,length(detections)]);

processFrames=p.processFrames;
if(isempty(processFrames))
    processFrames=1:frameNb;
end
frameNb=min([projData.frameNb,length(detections),length(processFrames)]);


%% create projection process saving independant projection location
% if(~isempty(p.process))
%   processRenderer = ProjRendering(processProj,p.name);
% end

colorIndx=p.colorIndx;
if(isempty(colorIndx))
    colorIndx=arrayfun(@(d) ones(1,length(d.zCoord(:,1))),detections,'unif',0);
end

colormap=p.colormap;
if(isempty(colormap))
    colormap=255*jet(length(unique(horzcat(colorIndx{:}))));
end

parfor fIdx=processFrames
    [XYProj,ZYProj,ZXProj,three]=processProj.loadFrame(1,fIdx);
    if(cumulative)
        detectionsAtFrame=detections;
        fColorIndx=p.colorIndx;
    else
        detectionsAtFrame=detections(fIdx);
        fColorIndx=colorIndx{fIdx};
    end
    % detectionsAtFrame.zCoord(:,1)=detectionsAtFrame.zCoord(:,1)/0.378;
    [overlayXY,overlayZY,overlayZX]=overlayProjDetections(XYProj,ZYProj,ZXProj, ...
        [projData.minXBorder projData.maxXBorder],[projData.minYBorder projData.maxYBorder],[projData.minZBorder projData.maxZBorder], ...
        detectionsAtFrame,colormap,fColorIndx,varargin{:});

    if(~isempty(p.process))
      p.process.saveFrame(1,fIdx,overlayXY,overlayZY,overlayZX);
    end
end

if(~isempty(p.process)) 
    ProjAnimation(p.process,'ortho').saveVideo([p.process.getOutputDir()  '.avi']);
end

