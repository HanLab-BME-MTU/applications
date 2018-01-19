function projImages=overlayProjTracksMovie(processProj,varargin)
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addRequired('processProj');
  ip.addOptional('tracksOrProcess',[]);
  ip.addOptional('tracks',[]);
  ip.addOptional('process',[]);
  ip.addOptional('processFrames',[]);
  ip.addOptional('dragonTail',[]);
  ip.addOptional('colormap',[]);
  ip.addOptional('colorIndx',[]);
  ip.addOptional('name','tracks');
  ip.parse(varargin{:});
  p=ip.Results;

  if(isempty(p.tracks))
    tracks=p.tracksOrProcess;
  else
    tracks=p.tracks;
  end
  
  if(isa(tracks,'Process'))
      tracksFinal=tracks.loadChannelOutput(1);
      tracks=TracksHandle(tracksFinal);
  end
      
  if(isempty(tracks))
      return;
  end
      
  projDataIdx=5;
  ref=[];
  try % Handle Project1D/ProjDyn with different implementation (current and then deprecated)
    ref=get(processProj,'ref');
    tracks=ref.applyBase(tracks,'');
    projData=processProj;
  catch
    try % Handle Project1D/ProjDyn different outFilePaths_spec (need to be defined through a class...)
      projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
    catch
      projDataIdx=4;
    end
    projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  end

  processFrames=p.processFrames;
  if(isempty(processFrames))
      processFrames=1:projData.frameNb;
  end

%% create projection process saving independant projection location
% if(~isempty(p.process))
% end

frameNb=projData.frameNb;
parfor fIdx=processFrames
  [XYProj,ZYProj,ZXProj,three]=processProj.loadFrame(1,fIdx);

  [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj, ...
      [projData.minXBorder projData.maxXBorder],[projData.minYBorder projData.maxYBorder],[projData.minZBorder projData.maxZBorder], ...
      fIdx,tracks,p.colormap,p.colorIndx,varargin{:});
  
  if(~isempty(p.process))
    p.process.saveFrame(1,fIdx,tracksXY,tracksZY,tracksZX);
  end
end

if(~isempty(p.process)) 
    ProjAnimation(p.process,'ortho').saveVideo([p.process.getOutputDir()  '.avi']);
end



