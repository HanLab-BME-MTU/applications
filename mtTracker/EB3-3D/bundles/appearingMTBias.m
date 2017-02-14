function appearingMTBias(kinTracks)


for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'appearingMTBias.');
    try
        kinTrack.addprop('distP1');
        kinTrack.addprop('distP2');
        kinTrack.addprop('MTP2Angle');
        kinTrack.addprop('MTP1Angle');
    catch
    end;
    kinTrack.MTP1Angle=[] ;
    kinTrack.distP1=[] ;
    kinTrack.MTP2Angle=[] ;
    kinTrack.distP2=[] ;   
    for mtIdx=1:length(kinTrack.appearingMTP1)
        mt=kinTrack.appearingMTP1(mtIdx);
        MTDisappFrame=mt.f(end);
        pIdx=find(kinTrack.f==MTDisappFrame);
        if(isempty(pIdx))
            pIdx=length(kinTrack.f);
        end
        %% Angle to kinPole axis
        mtVectorP1KinRef=[mt.pole1.x(end); mt.pole1.y(end); mt.pole1.z(end) ];
        kinVectorP1KinRef=[kinTrack.pole1.x(pIdx); kinTrack.pole1.y(pIdx); kinTrack.pole1.z(pIdx)];
        MTAnglesP1Kin=vectorAngleND(mtVectorP1KinRef,kinVectorP1KinRef);
        distP1=sin(MTAnglesP1Kin)*mt.rho(end);
        kinTrack.MTP1Angle=[kinTrack.MTP1Angle; MTAnglesP1Kin] ;
        kinTrack.distP1=[kinTrack.distP1; distP1] ;
    end
    for mtIdx=1:length(kinTrack.appearingMTP2)
        mt=kinTrack.appearingMTP2(mtIdx);
        MTDisappFrame=mt.f(end);
        pIdx=find(kinTrack.f==MTDisappFrame);
        if(isempty(pIdx))
            pIdx=length(kinTrack.f);
        end
        %% Angle to kinPole axis
        mtVectorP2KinRef=[mt.pole2.x(end); mt.pole2.y(end); mt.pole2.z(end) ];
        kinVectorP2KinRef=[kinTrack.pole2.x(pIdx); kinTrack.pole2.y(pIdx); kinTrack.pole2.z(pIdx)];
        MTAnglesP2Kin=vectorAngleND(mtVectorP2KinRef,kinVectorP2KinRef);
        kinTrack.MTP2Angle=[kinTrack.MTP2Angle; MTAnglesP2Kin] ;
        distP2=sin(MTAnglesP2Kin)*mt.rho(end);
        kinTrack.distP2=[kinTrack.distP2; distP2] ;        
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