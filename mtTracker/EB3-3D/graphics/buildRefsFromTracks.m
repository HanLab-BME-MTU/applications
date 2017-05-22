function refs=buildRefsFromTracks(originsTracksOrProcess,ZTracksOrProcess,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('originsTracksOrProcess');
ip.addRequired('ZTracksOrProcess');
ip.addParamValue('process', []);
ip.parse(originsTracksOrProcess,ZTracksOrProcess,varargin{:});
p=ip.Results;

if(isa(originsTracksOrProcess,'Process'))
    tmp=load(originsTracksOrProcess.outFilePaths_{1}); origins=tmp.tracks;
else
    origins=originsTracksOrProcess;
end

if(isa(ZTracksOrProcess,'Process'))
    tmp=load(ZTracksOrProcess.outFilePaths_{1}); ZTracks=tmp.tracks;
else
    ZTracks=ZTracksOrProcess;
end

refs(length(origins),length(ZTracks))=FrameOfRef();
for zIdx=1:length(ZTracks)
    for orIdx=1:length(origins)
        ref=FrameOfRef();
        ref.setOriginFromTrack(origins(orIdx));
        ref.setZFromTrack(ZTracks(zIdx));
        ref.genBaseFromZ();
        refs(orIdx,zIdx)=ref;
    end
end

process=ip.Results.process;
if(~isempty(process))
    outputDir=[process.getOwner().outputDirectory_ filesep 'refs'];
    mkdirRobust(outputDir);
    save([outputDir filesep 'refs.mat'],'refs');
    process.setOutFilePaths({[outputDir filesep 'refs.mat']})
    pa = process.getParameters();
    pa.parameters = ip.Results;
    process.setParameters(pa);
    process.setDateTime();
end

