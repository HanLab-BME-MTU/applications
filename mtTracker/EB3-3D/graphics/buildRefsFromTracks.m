function refs=buildRefsFromTracks(originsTracksOrProcess,ZTracksOrProcess)
if(isa(originsTracksOrProcess,'Process'))
    tmp=load(originsTracksOrProcess.outFilePaths_{1}); origins=tmp.fiducialsTracks;
else
    origins=originsTracksOrProcess;
end

if(isa(ZTracksOrProcess,'Process'))
    tmp=load(ZTracksOrProcess.outFilePaths_{1}); ZTracks=tmp.fiducialsTracks;
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
