function [kinTracksP1,kinTracksP2]=mapMTApparitionToKin(kinTracksP1,kinTracksP2,EB3TracksP1,EB3TracksP2,angleCutoff)
  
    EB3StartFrames=[EB3TracksP1.startFrame]';
    EB3RhoP1=arrayfun(@(t) t.rho(1),EB3TracksP1);
    EB3RhoP2=arrayfun(@(t) t.rho(1),EB3TracksP2);

    for kinIdx=1:length(kinTracksP1)
        kinTrackP1=kinTracksP1(kinIdx);
        kinTrackP2=kinTracksP2(kinIdx);
        progressText(kinIdx/length(kinTracksP1),'mapMTApparitionToKin.');
        try
            kinTrackP1.addprop('appearingMT');
            kinTrackP2.addprop('appearingMT');            
        catch
        end;
        kinTrackP1.appearingMT=[];
        kinTrackP2.appearingMT=[];       
        for pIdx=1:length(kinTrackP1.f)
            fIdx=kinTrackP1.f(pIdx);
            coexistingEB3=(EB3StartFrames==fIdx); % EB3 and Kin Co-exist if EB3 appear when Kin is alive.
            if(any(coexistingEB3))
                % associate EB3 to kin
                % Compute angle 
                % threshold manually
                
                % Kinetochore distance from each pole at frame <fIdx>
                kinRhoP1=kinTrackP1.rho(pIdx);
                kinRhoP2=kinTrackP2.rho(pIdx);

                % Associate MT appearance to kinPole refererial using
                % Distance to pole
                P1KinAssociatedMT=(EB3RhoP1<kinRhoP1)&coexistingEB3;
                P2KinAssociatedMT=(EB3RhoP2<kinRhoP2)&coexistingEB3;
                P1KinAssociatedMTIndex=find(P1KinAssociatedMT);
                P2KinAssociatedMTIndex=find(P2KinAssociatedMT);
                
                %% Angle to kinPole axis          
                MTVectorP1KinRef=cell2mat(arrayfun(@(t) [t.x(1);t.y(1);t.z(1)],EB3TracksP1(P1KinAssociatedMT),'unif',0)');
                MTVectorP2KinRef=cell2mat(arrayfun(@(t) [t.x(1);t.y(1);t.z(1)],EB3TracksP2(P2KinAssociatedMT),'unif',0)');

                kinVectorP1KinRef=[kinTrackP1.x(pIdx); kinTrackP1.y(pIdx); kinTrackP1.z(pIdx) ];
                kinVectorP2KinRef=[kinTrackP2.x(pIdx); kinTrackP2.y(pIdx); kinTrackP2.z(pIdx) ];
                %%                
                MTAnglesP1Kin= vectorAngleND(MTVectorP1KinRef,kinVectorP1KinRef);
                MTAnglesP2Kin= vectorAngleND(MTVectorP2KinRef,kinVectorP2KinRef);
                %%
                appearingMTP1Kin=EB3TracksP1(P1KinAssociatedMTIndex(abs(MTAnglesP1Kin)<angleCutoff));
                appearingMTP2Kin=EB3TracksP2(P2KinAssociatedMTIndex(abs(MTAnglesP2Kin)<angleCutoff));
                              
                if(~isempty(appearingMTP1Kin))
                    kinTrackP1.appearingMT=[kinTrackP1.appearingMT; appearingMTP1Kin] ;
                end
                if(~isempty(appearingMTP2Kin))
                    kinTrackP2.appearingMT=[kinTrackP2.appearingMT; appearingMTP2Kin] ;
                end                
            end
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