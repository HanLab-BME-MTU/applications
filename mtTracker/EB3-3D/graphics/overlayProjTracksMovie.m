function projImages=overlayProjTracksMovie(processSingleProj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('processSingleProj');
ip.addOptional('tracks',[]);
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
system(['mkdir -p ' fileparts(savePath)]);
frameNb=projData.frameNb;
for fIdx=1:frameNb
    XYProj=imread(sprintf(XYProjTemplate,fIdx));
    ZYProj=imread(sprintf(ZYProjTemplate,fIdx));
    ZXProj=imread(sprintf(ZXProjTemplate,fIdx));
    [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj, ...
      [projData.minXBorder projData.maxXBorder],[projData.minYBorder projData.maxYBorder],[projData.minZBorder projData.maxZBorder], ...
      fIdx,tracks,p.colormap,p.colorIndx);

    %% Use Z to index image line (going up)
%     tracksXY=permute(tracksXY,[2 1 3]);
    tracksZY=permute(tracksZY,[2 1 3]);
    tracksZX=permute(tracksZX,[2 1 3]);
    three=projMontage(tracksXY,tracksZX,tracksZY);
    imwrite(three,sprintf(savePath,fIdx));
end

% save as video 
video = VideoWriter([fileparts(processSingleProj.outFilePaths_{4})  filesep  p.name '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75
open(video)
for fIdx=1:frameNb
  three=[imread(sprintf(savePath,fIdx))];
    writeVideo(video,three);
end
close(video)
