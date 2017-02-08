function [kinTracks]=mapMTApparitionToKin(kinTracks,EB3Tracks,angleCutoff)
  
    EB3StartFrames=[EB3Tracks.startFrame]';
    EB3RhoP1=arrayfun(@(t) t.poleRef{1}.rho(1),EB3Tracks);
    EB3RhoP2=arrayfun(@(t) t.poleRef{2}.rho(1),EB3Tracks);

    for kinIdx=1:length(kinTracks)
        kinTrack=kinTracks(kinIdx);
        progressText(kinIdx/length(kinTracks),'mapMTApparitionToKin.');
        try
            kinTrack.addprop('appearingMT');
        catch
        end;
%        kinTrack.appearingMT=cell(1,length(kinTrack.poleRef));
        kinTrack.appearingMT=[];

        for pIdx=1:length(kinTrack.f)
            fIdx=kinTrack.f(pIdx);
            coexistingEB3=(EB3StartFrames==fIdx); % EB3 and Kin Co-exist if EB3 appear when Kin is alive.
            if(any(coexistingEB3))
                % associate EB3 to kin
                % Compute angle 
                % threshold manually
                
                % Kinetochore distance from each pole at frame <fIdx>
                kinRhoP1=kinTrack.poleRef{1}.rho(pIdx);
                kinRhoP2=kinTrack.poleRef{2}.rho(pIdx);

                % Associate MT appearance to kinPole refererial using
                % Distance to pole
                P1KinAssociatedMT=(EB3RhoP1<kinRhoP1)&coexistingEB3;
                P2KinAssociatedMT=(EB3RhoP2<kinRhoP2)&coexistingEB3;
                P1KinAssociatedMTIndex=find(P1KinAssociatedMT);
                P2KinAssociatedMTIndex=find(P2KinAssociatedMT);
                
                %% Angle to kinPole axis          
                MTVectorP1KinRef=cell2mat(arrayfun(@(t) [t.poleRef{1}.x(1);t.poleRef{1}.y(1);t.poleRef{1}.z(1)],EB3Tracks(P1KinAssociatedMT),'unif',0)');
                MTVectorP2KinRef=cell2mat(arrayfun(@(t) [t.poleRef{2}.x(1);t.poleRef{2}.y(1);t.poleRef{2}.z(1)],EB3Tracks(P2KinAssociatedMT),'unif',0)');

                kinVectorP1KinRef=[kinTrack.poleRef{1}.x(pIdx); kinTrack.poleRef{1}.y(pIdx); kinTrack.poleRef{1}.z(pIdx) ];
                kinVectorP2KinRef=[kinTrack.poleRef{2}.x(pIdx); kinTrack.poleRef{2}.y(pIdx); kinTrack.poleRef{2}.z(pIdx) ];
                
                %%                
                MTAnglesP1Kin= vectorAngleND(MTVectorP1KinRef,kinVectorP1KinRef);
                MTAnglesP2Kin= vectorAngleND(MTVectorP2KinRef,kinVectorP2KinRef);
                %%
                appearingMTP1Kin=EB3Tracks(P1KinAssociatedMTIndex(abs(MTAnglesP1Kin)<angleCutoff));
                appearingMTP2Kin=EB3Tracks(P2KinAssociatedMTIndex(abs(MTAnglesP2Kin)<angleCutoff));
                              
                if(~isempty(appearingMTP1Kin))
                    kinTrack.appearingMT=[kinTrack.appearingMT; appearingMTP1Kin] ;
                end
                if(~isempty(appearingMTP2Kin))
                    kinTrack.appearingMT=[kinTrack.appearingMT; appearingMTP2Kin] ;
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