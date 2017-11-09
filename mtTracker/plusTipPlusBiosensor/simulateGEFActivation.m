function renderGEFDiffSimulation(MD,tracks,diffCoef,diffTime,varargin)%,,diffCoeff,varargin)
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParamValue('showAll',false,@islogical);
ip.parse( varargin{:});
p=ip.Results;

imageDiff=zeros(1,MD.nFrames_-1);
parfor i=1:(MD.nFrames_-1)
    imageDiff(i)=sum(sum(abs(MD.getChannel(1).loadImage(i+1)-MD.getChannel(1).loadImage(i))));
end

jumpTresh=median(imageDiff)+3*1.4826*mad(imageDiff,1);
jump=find(imageDiff>jumpTresh)+1;

if(p.showAll)
    plot(imageDiff)
    vline(jump,'r-');
end

endOnJumpIdx=ismember([tracks.endFrame],jump);
uncutTracks=tracks(endOnJumpIdx);
amiraWriteTracks([MD.outputDirectory_ filesep 'activationSimulation' filesep 'uncutTracks/uncutTracks.am'],uncutTracks);

activationMap=cell(1,MD.nFrames_);
activationMap=cellfun(@(x) zeros(MD.imSize_),activationMap,'unif',false);

for i=1:length(uncutTracks)
    tr=uncutTracks(i);
    for t=tr.endFrame:min(MD.nFrames_,(tr.endFrame+diffTime))
        idx=t-tr.endFrame+1;
        activationMap{t}=MidpointCircle(activationMap{t},diffCoef*idx+1,tr.y(end), tr.x(end),diffTime-idx);
    end
end
%%
mkdir([MD.outputDirectory_ filesep 'activationSimulation' filesep 'activationMap']);
parfor i=1:(MD.nFrames_)
   imwrite(uint8(activationMap{i}),[MD.outputDirectory_ filesep 'activationSimulation' filesep 'activationMap/map' num2str(i,'%04d') '.tif']);
end