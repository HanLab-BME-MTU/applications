function projImages=overlayProjDetectionMovie(processSingleProj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('processSingleProj');
ip.addOptional('detections',[]);
ip.addOptional('colormap',[]);
ip.addOptional('process',[]);
ip.addOptional('colorIndx',[]);
ip.addOptional('name','detections');
ip.parse(varargin{:});
p=ip.Results;

detections=p.detections;

%% testing imwarp to crop the image
XYProjTemplate=processSingleProj.outFilePaths_{1};
ZYProjTemplate=processSingleProj.outFilePaths_{2};
ZXProjTemplate=processSingleProj.outFilePaths_{3};
projData=load(processSingleProj.outFilePaths_{4},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
savePath=[fileparts(processSingleProj.outFilePaths_{4}) filesep p.name filesep 'frame_nb%04d.png'];
mkdirRobust([ fileparts(savePath)]);
frameNb=projData.frameNb;

%% create projection process saving independant projection location

if(~isempty(p.process)) 
    p.process.setOutFilePaths({savePath});
end




for fIdx=1:frameNb
    XYProj=imread(sprintfPath(XYProjTemplate,fIdx));
    ZYProj=imread(sprintfPath(ZYProjTemplate,fIdx));
    ZXProj=imread(sprintfPath(ZXProjTemplate,fIdx));
    detectionsAtFrame=detections(fIdx);
    [tracksXY,tracksZY,tracksZX]=overlayProjDetections(XYProj,ZYProj,ZXProj, ...
      [projData.minXBorder projData.maxXBorder],[projData.minYBorder projData.maxYBorder],[projData.minZBorder projData.maxZBorder], ...
      detectionsAtFrame,p.colormap,p.colorIndx{fIdx},varargin{:});

    %% Use Z to index image line (going up)
%     tracksXY=permute(tracksXY,[2 1 3]);
    tracksZY=permute(tracksZY,[2 1 3]);
    tracksZX=permute(tracksZX,[2 1 3]);
    three=projMontage(tracksXY,tracksZX,tracksZY);
    imwrite(three,sprintfPath(savePath,fIdx));
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

