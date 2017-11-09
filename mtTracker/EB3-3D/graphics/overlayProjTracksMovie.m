function projImages=overlayProjTracksMovie(processSingleProj,varargin)
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addRequired('processSingleProj');
  ip.addOptional('tracksOrProcess',[]);
  ip.addOptional('tracks',[]);
  ip.addOptional('process',[]);
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
      
%% testing imwarp to crop the image
XYProjTemplate=processSingleProj.outFilePaths_{1};
ZYProjTemplate=processSingleProj.outFilePaths_{2};
ZXProjTemplate=processSingleProj.outFilePaths_{3};
projDataIdx=5;
try % Handle Project1D/ProjDyn different outFilePaths_spec (need to be defined through a class...)
  projData=load(processSingleProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
catch
  projDataIdx=4;
end
projData=load(processSingleProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');

savePath=[fileparts(processSingleProj.outFilePaths_{projDataIdx}) filesep p.name filesep 'frame_nb%04d.png'];
outputDir=fileparts(savePath);
mkdirRobust(outputDir);

%% create projection process saving independant projection location
if(~isempty(p.process))
  mkdirRobust([outputDir filesep 'XY'])
  mkdirRobust([outputDir filesep 'YZ'])
  mkdirRobust([outputDir filesep 'XZ']) 
  mkdirRobust([outputDir filesep 'three']) 
  load(processSingleProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  save([outputDir filesep 'limits.mat'],'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  p.process.setOutFilePaths({[outputDir filesep 'XY' filesep 'frame_nb%04d.png'], ...
    [outputDir filesep 'YZ' filesep 'frame_nb%04d.png'], ...
    [outputDir filesep 'XZ' filesep 'frame_nb%04d.png'],...
    [outputDir filesep 'three' filesep 'frame_nb%04d.png'],...
    [outputDir filesep 'limits.mat']});
end

frameNb=projData.frameNb;
parfor fIdx=1:frameNb
  XYProj=imread(sprintfPath(XYProjTemplate,fIdx));
  ZYProj=imread(sprintfPath(ZYProjTemplate,fIdx));
  ZXProj=imread(sprintfPath(ZXProjTemplate,fIdx));
  [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj, ...
      [projData.minXBorder projData.maxXBorder],[projData.minYBorder projData.maxYBorder],[projData.minZBorder projData.maxZBorder], ...
      fIdx,tracks,p.colormap,p.colorIndx,varargin{:});
  
  %% Use Z to index image line (going up)
  %     tracksXY=permute(tracksXY,[2 1 3]);
  if(~isempty(p.process))
      % save the maximum intensity projections
      imwrite(tracksXY, sprintfPath(p.process.outFilePaths_{1},fIdx), 'Compression', 'none');
      imwrite(tracksZY, sprintfPath(p.process.outFilePaths_{2},fIdx), 'Compression', 'none');
      imwrite(tracksZX, sprintfPath(p.process.outFilePaths_{3},fIdx), 'Compression', 'none');
  end
  tracksZY=permute(tracksZY,[2 1 3]);
  tracksZX=permute(tracksZX,[2 1 3]);
  three=projMontage(tracksXY,tracksZX,tracksZY);
  imwrite(three,sprintfPath(p.process.outFilePaths_{4},fIdx));
end

% save as video
video = VideoWriter([fileparts(processSingleProj.outFilePaths_{projDataIdx})  '-'  p.name '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75
open(video)
for fIdx=1:frameNb
  three=[imread(sprintfPath(p.process.outFilePaths_{4},fIdx))];
    writeVideo(video,three);
end
close(video)

