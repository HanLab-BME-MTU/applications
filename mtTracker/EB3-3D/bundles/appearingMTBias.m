function appearingMTBias(kinTracks)

for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'appearingMTBias.');
    try
        kinTrack.addprop('MTP2Angle');
        kinTrack.addprop('MTP1Angle');        
    catch
    end;
    for mtIdx=1:length(kinTrack.appearingMTP1)
        mt=kinTrack.appearingMTP1(mtIdx);
        MTDisappFrame=mt.f(end);
        pIdx=find(kinTrack.f==MTDisappFrame);
        if(isempty(pIdx))
            pIdx=length(kinTrack.f);
        end
        %% Angle to kinPole axis
        mtVectorP1KinRef=[mt.poleRef{1}.x(end); mt.poleRef{1}.y(end); mt.poleRef{1}.z(end) ];
        kinVectorP1KinRef=[kinTrack.poleRef{1}.x(pIdx); kinTrack.poleRef{1}.y(pIdx); kinTrack.poleRef{1}.z(pIdx)];
        MTAnglesP1Kin=vectorAngleND(mtVectorP1KinRef,kinVectorP1KinRef);
        kinTrack.MTP1Angle=[kinTrack.MTP1Angle; MTAnglesP1Kin] ;
    end
    for mtIdx=1:length(kinTrack.appearingMTP2)
        mt=kinTrack.appearingMTP2(mtIdx);
        MTDisappFrame=mt.f(end);
        pIdx=find(kinTrack.f==MTDisappFrame);
        if(isempty(pIdx))
            pIdx=length(kinTrack.f);
        end
        %% Angle to kinPole axis
        mtVectorP2KinRef=[mt.poleRef{2}.x(end); mt.poleRef{2}.y(end); mt.poleRef{2}.z(end) ];
        kinVectorP2KinRef=[kinTrack.poleRef{2}.x(pIdx); kinTrack.poleRef{2}.y(pIdx); kinTrack.poleRef{2}.z(pIdx)];
        MTAnglesP2Kin=vectorAngleND(mtVectorP2KinRef,kinVectorP2KinRef);
        kinTrack.MTP2Angle=[kinTrack.MTP2Angle; MTAnglesP2Kin] ;
    end    
end
end

function angle=vectorAngleND(a,b)
    if (size(a,2)==1)
        a=repmat(a,1,size(b,2));
    end
    if (size(b,2)==1)
        b=repmat(b,1,size(a,2));
    end
    f=@(a,b)atan2(norm(cross(a,b)), dot(a,b));
    angle=arrayfun(@(i) f(a(:,i),b(:,i)),1:size(a,2));
end