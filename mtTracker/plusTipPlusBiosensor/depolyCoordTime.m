function [Pr,XY,time]=depolyCoordTime(MD,varargin)%,,diffCoeff,varargin)
ip=inputParser(); ip.CaseSensitive = false; ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addOptional('tracks',[],@(x) (isstruct(x)));
ip.addParamValue('showAll',false,@islogical);
ip.addParamValue('chIdx',1,@isnumeric);
ip.parse(MD,varargin{:});
p=ip.Results;

outputFolder=[MD.outputDirectory_ filesep 'MT-GEF' filesep 'MTDepoly' filesep]; mkdir(outputFolder);
outputDepolyCoordFile=[outputFolder filesep 'XYTCoordinate.mat'];
outputTrackFile=[outputFolder filesep 'uninterruptedTracks.mat'];
outputTrackAmiraFiles=[outputFolder filesep 'uninterruptedTracksAmira' filesep 'uninterruptedTracks.am'];

Pr=ExternalProcess(MD,'depolyCoordTime');
Pr.setInFilePaths({});
Pr.setOutFilePaths({outputDepolyCoordFile,outputTrackFile,outputTrackAmiraFiles});
pa = Pr.getParameters();
pa.parameters = p;
Pr.setParameters(pa);
Pr.setDateTime();
MD.addProcess(Pr);

if(isempty(p.tracks))
    trackingProcessIdx=MD.getProcessIndex('TrackingProcess',1,1);
    if isempty(trackingProcessIdx)
        error('Tracking has not been run! Please run tracking or input a valid TrackHandle object as a second argument.');
    end
    tracksFinal=MD.getProcess(trackingProcessIdx).loadChannelOutput(1);
    tracks=TracksHandle(tracksFinal);
else
    tracks=p.tracks;
end

imageDiff=zeros(1,MD.nFrames_-1);
parfor i=1:(MD.nFrames_-1)
    imageDiff(i)=sum(sum(abs(MD.getChannel(1).loadImage(i+1)-MD.getChannel(1).loadImage(i))));
end

jumpTresh=median(imageDiff)+3*1.4826*mad(imageDiff,1);
beforeJumpFrame=find(imageDiff>jumpTresh);

%%
if(p.showAll)
  %%
    plot(imageDiff)
    vline(beforeJumpFrame,'r-');
end
%%
endOnJumpIdx=ismember([tracks.endFrame],beforeJumpFrame);
uninterruptedTracks=tracks(~endOnJumpIdx);
save(outputTrackFile,'uninterruptedTracks');
amiraWriteTracks(outputTrackAmiraFiles,uninterruptedTracks);


XY=zeros(length(uninterruptedTracks),2);
time=zeros(length(uninterruptedTracks),1);
for i=1:length(uninterruptedTracks)
    tr=uninterruptedTracks(i);
    time(i)=tr.t(end);
    XY(i,:)=[tr.x(end) tr.y(end)];
end
save(outputDepolyCoordFile,'XY','time');