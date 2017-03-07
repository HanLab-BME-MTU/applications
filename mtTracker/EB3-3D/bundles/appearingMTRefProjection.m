function appearingMTRefProjection(kinTracks,MTFieldNameP1,MTFieldNameP2)


for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'appearingMTBias.');
    newRefFieldP1=[MTFieldNameP1 'KinRef'];
    newRefFieldP2=[MTFieldNameP2 'KinRef'];

    try
        kinTrack.addprop(newRefFieldP1);
        kinTrack.addprop(newRefFieldP2);
    catch
    end;
    setfield(kinTrack,newRefFieldP1,[]);
    setfield(kinTrack,newRefFieldP2,[]);
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
    
    P1AssoMT=getfield(kinTrack,MTFieldNameP1);
    appearingMTKinP1Ref=[];
    for mtIdx=1:length(P1AssoMT)
        mt=P1AssoMT(mtIdx);
        appearingMTKinP1Ref=[appearingMTKinP1Ref; refP1.applyBase(mt,'')];
    end
    setfield(kinTrack,newRefFieldP1,appearingMTKinP1Ref);
    
    P2AssoMT=getfield(kinTrack,MTFieldNameP2);
    appearingMTKinP2Ref=[];
    for mtIdx=1:length(P2AssoMT)
        mt=P2AssoMT(mtIdx);
        appearingMTKinP2Ref=[appearingMTKinP2Ref; refP2.applyBase(mt,'')];
    end
    setfield(kinTrack,newRefFieldP2,appearingMTKinP2Ref);

end
end

