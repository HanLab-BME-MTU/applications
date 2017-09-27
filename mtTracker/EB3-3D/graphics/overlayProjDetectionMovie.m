function overlayProjDetectionMovie(processSingleProj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('processSingleProj');
ip.addOptional('detections',[]);
ip.addOptional('colormap',[]);
ip.addOptional('process',[]);
ip.addOptional('processFrames',[]);
ip.addOptional('cumulative',false);
ip.addOptional('colorIndx',[]);
ip.addOptional('name','detections');
ip.parse(varargin{:});
p=ip.Results;

detections=p.detections;
cumulative=p.cumulative;
%% testing imwarp to crop the image
XYProjTemplate=processSingleProj.outFilePaths_{1};
ZYProjTemplate=processSingleProj.outFilePaths_{2};
ZXProjTemplate=processSingleProj.outFilePaths_{3};
projData=load(processSingleProj.outFilePaths_{5},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
savePath=[fileparts(processSingleProj.outFilePaths_{5}) filesep p.name filesep 'frame_nb%04d.png'];
outputDir=fileparts(savePath);
mkdirRobust([fileparts(savePath)]);
frameNb=min([projData.frameNb,length(detections)]);

processFrames=p.processFrames;
if(isempty(processFrames))
    processFrames=1:frameNb;
end
frameNb=min([projData.frameNb,length(detections),length(processFrames)]);

%% create projection process saving independant projection location
if(~isempty(p.process))
  mkClrDir([outputDir filesep 'XY'])
  mkClrDir([outputDir filesep 'YZ'])
  mkClrDir([outputDir filesep 'XZ']) 
  mkClrDir([outputDir filesep 'three'])  
  minXBorder=projData.minXBorder;maxXBorder=projData.maxXBorder;
  minYBorder=projData.minYBorder;maxYBorder=projData.maxYBorder;
  minZBorder=projData.minZBorder;maxZBorder=projData.maxZBorder;
  save([outputDir filesep 'limits.mat'],'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  p.process.setOutFilePaths({[outputDir filesep 'XY' filesep 'frame_nb%04d.png'], ...
    [outputDir filesep 'YZ' filesep 'frame_nb%04d.png'], ...
    [outputDir filesep 'XZ' filesep 'frame_nb%04d.png'],...
    [outputDir filesep 'three' filesep 'frame_nb%04d.png'],...
    [outputDir filesep 'limits.mat']});
end

if(p.cumulative)
    tmpdetections=copy(detections(1));
    for i=2:length(detections)
        tmpdetections.concatenate(detections(i));
    end
    detections=tmpdetections;
end
    

parfor fIdx=processFrames
    XYProj=imread(sprintfPath(XYProjTemplate,fIdx));
    ZYProj=imread(sprintfPath(ZYProjTemplate,fIdx));
    ZXProj=imread(sprintfPath(ZXProjTemplate,fIdx));
    if(cumulative)
        detectionsAtFrame=detections;
    else
        detectionsAtFrame=detections(fIdx);
    end
    [tracksXY,tracksZY,tracksZX]=overlayProjDetections(XYProj,ZYProj,ZXProj, ...
        [projData.minXBorder projData.maxXBorder],[projData.minYBorder projData.maxYBorder],[projData.minZBorder projData.maxZBorder], ...
        detectionsAtFrame,p.colormap,p.colorIndx{fIdx},varargin{:});
    
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
    if(~isempty(p.process))
        imwrite(three, sprintfPath(p.process.outFilePaths_{4},fIdx), 'Compression', 'none');
    end    
end



% save as video
video = VideoWriter([fileparts(processSingleProj.outFilePaths_{4})  '-'  p.name '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75
open(video)
for fIdx=1:frameNb
  three=[imread(sprintfPath(p.process.outFilePaths_{4},fIdx))];
    writeVideo(video,three);
end
close(video)

