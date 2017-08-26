function projImages=overlayProjTracksMovie(processSingleProj,varargin)
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addRequired('processSingleProj');
  ip.addOptional('tracks',[]);
  ip.addOptional('process',[]);
  ip.addOptional('colormap',[]);
  ip.addOptional('colorIndx',[]);
  ip.addOptional('name','tracks');
  ip.parse(varargin{:});
  p=ip.Results;

  tracks=p.tracks;


%% testing imwarp to crop the image
XYProjTemplate=processSingleProj.outFilePaths_{1};
ZYProjTemplate=processSingleProj.outFilePaths_{2};
ZXProjTemplate=processSingleProj.outFilePaths_{3};
projData=load(processSingleProj.outFilePaths_{4},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
savePath=[fileparts(processSingleProj.outFilePaths_{4}) filesep p.name filesep 'frame_nb%04d.png'];
outputDir=fileparts(savePath);
mkdirRobust(outputDir);

%% create projection process saving independant projection location
if(~isempty(p.process))
  mkdirRobust([outputDir filesep 'XY'])
  mkdirRobust([outputDir filesep 'YZ'])
  mkdirRobust([outputDir filesep 'XZ'])    
  p.process.setOutFilePaths({[outputDir filesep 'XY' filesep 'frame_nb%04d.png'], ...
    [outputDir filesep 'YZ' filesep 'frame_nb%04d.png'], ...
    [outputDir filesep 'XZ' filesep 'frame_nb%04d.png'],...
    [outputDir filesep 'limits.mat']});
    load(processSingleProj.outFilePaths_{4},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
    save([outputDir filesep 'limits.mat'],'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
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
  tracksZY=permute(tracksZY,[2 1 3]);
  tracksZX=permute(tracksZX,[2 1 3]);
  three=projMontage(tracksXY,tracksZX,tracksZY);
  imwrite(three,sprintfPath(savePath,fIdx));
  if(~isempty(p.process))
    % save the maximum intensity projections
    imwrite(tracksXY, sprintfPath(p.process.outFilePaths_{1},fIdx), 'Compression', 'none');
    imwrite(tracksZY, sprintfPath(p.process.outFilePaths_{2},fIdx), 'Compression', 'none');
    imwrite(tracksZX, sprintfPath(p.process.outFilePaths_{3},fIdx), 'Compression', 'none');
  end
end

% save as video
video = VideoWriter([fileparts(processSingleProj.outFilePaths_{4})  '-'  p.name '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75
open(video)
for fIdx=1:frameNb
  three=[imread(sprintfPath(savePath,fIdx))];
    writeVideo(video,three);
end
close(video)

