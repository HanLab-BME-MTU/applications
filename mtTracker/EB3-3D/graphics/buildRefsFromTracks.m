function refs=buildRefsFromTracks(origins,ZTracks)
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
