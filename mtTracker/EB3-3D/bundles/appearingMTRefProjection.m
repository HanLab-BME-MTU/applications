function appearingMTRefProjection(kinTracks)


for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'appearingMTBias.');
    try
        kinTrack.addprop('appearingMTKinP1Ref');
        kinTrack.addprop('appearingMTKinP2Ref');
    catch
    end;
    kinTrack.appearingMTKinP1Ref=[] ;
    kinTrack.appearingMTKinP2Ref=[] ;
    refP1=FrameOfRef();
    refP1.origin=(kinTrack.pole1.ref.origin(kinTrack.f,:));
    refP1.setZFromTrack(kinTrack);
    refP1.genBaseFromZ();
    refP1.applyBaseToTrack(kinTrack,'KP1');

    refP2=FrameOfRef();
    refP2.origin=(kinTrack.pole2.ref.origin(kinTrack.f,:));
    refP2.setZFromTrack(kinTrack);
    refP2.genBaseFromZ();
    refP2.applyBaseToTrack(kinTrack,'KP2');

    for mtIdx=1:length(kinTrack.appearingMTP1)
        mt=kinTrack.appearingMTP1(mtIdx);
        kinTrack.appearingMTKinP1Ref=[kinTrack.appearingMTKinP1Ref; refP1.applyBaseToTrack(mt,'')];
    end
    for mtIdx=1:length(kinTrack.appearingMTP2)
        mt=kinTrack.appearingMTP2(mtIdx);
        kinTrack.appearingMTKinP2Ref=[kinTrack.appearingMTKinP2Ref; refP2.applyBaseToTrack(mt,'')];
    end
end
end

