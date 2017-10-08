function [poleRefs,P1,P2]=buildSpindleRef(poleMovieInfo,pixelSize,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('poleMovieInfo',[],@(s) isstruct(s));
ip.addOptional('pixelSize',[]);
ip.addOptional('processDetectPoles',[]);
ip.addParameter('process',[]);
ip.parse(poleMovieInfo,pixelSize,varargin{:});
p=ip.Results;
% WARNING: this is not a trajectory, merely a collection of poles to ease
% implementation.
% pixelSize=p.pixelSize;
if(~isempty(p.processDetectPoles))
    tmp=load(p.processDetectPoles.outFilePaths_{1});
    poleMovieInfo=tmp.poleMovieInfo;
    MD=p.processDetectPoles.getOwner();
    pixelSize=MD.pixelSize_;
end

P1=TracksHandle();
P1.x=arrayfun(@(d) pixelSize*(d.xCoord(1,1)-1)+1,poleMovieInfo)';
P1.y=arrayfun(@(d) pixelSize*(d.yCoord(1,1)-1)+1,poleMovieInfo)';
P1.z=arrayfun(@(d) pixelSize*(d.zCoord(1,1)-1)+1,poleMovieInfo)';
P1.endFrame=length(poleMovieInfo);
P1.startFrame=1;

P2=TracksHandle();
P2.x=arrayfun(@(d) pixelSize*(d.xCoord(2,1)-1)+1,poleMovieInfo)';
P2.y=arrayfun(@(d) pixelSize*(d.yCoord(2,1)-1)+1,poleMovieInfo)';
P2.z=arrayfun(@(d) pixelSize*(d.zCoord(2,1)-1)+1,poleMovieInfo)';
P2.endFrame=length(poleMovieInfo);
P2.startFrame=1;

refP1=FrameOfRef();
refP1.setOriginFromTrack(P1);
refP1.setZFromTrack(P2);
refP1.genBaseFromZ();

refP2=FrameOfRef();
refP2.setOriginFromTrack(P2);
refP2.setZFromTrack(P1);
refP2.genBaseFromZ();

poleRefs=[refP1 refP2];

process=p.process;
if(~isempty(process))
    outputDirPoleDetect=[fileparts(p.processDetectPoles.outFilePaths_{1}) filesep];
    mkdir(outputDirPoleDetect);
    save([outputDirPoleDetect filesep 'refs.mat'],'poleRefs');
    process.setOutFilePaths({[outputDirPoleDetect filesep 'refs.mat']});
    pa = process.getParameters();
    pa.parameters = ip.Results;
    process.setParameters(pa);
    process.setDateTime();
end
