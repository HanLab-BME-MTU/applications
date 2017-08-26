function mapMTTipsToKin(kinTracks,EB3Tracks,cutoff,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('distType','normalDist');
ip.addParameter('minDistFromPole',500);
ip.parse(varargin{:});
p=ip.Results;

EB3StartFrames=[EB3Tracks.startFrame]';
EB3EndFrames=[EB3Tracks.endFrame]';

% Pole info to circumvent projecion bug
P1P2Ref=kinTracks(1).pole1.ref;
P2P1Ref=kinTracks(1).pole2.ref;

for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'mapMTApparitionToKin.');
    try
        kinTrack.addprop('associatedTipsP2Idx');
        kinTrack.addprop('associatedTipsP1Idx');
        kinTrack.addprop('associatedTipsP1');
        kinTrack.addprop('associatedTipsP2');
        kinTrack.addprop('distTipsP1');
        kinTrack.addprop('distTipsP2');
    catch
    end;
    %        kinTrack.appearingMT=cell(1,length(kinTrack.poleRef));
    kinTrack.associatedTipsP1=[];
    kinTrack.associatedTipsP2=[];
    kinTrack.associatedTipsP1Idx=[];
    kinTrack.associatedTipsP2Idx=[];    
    kinTrack.distTipsP1=[];
    kinTrack.distTipsP2=[];
    for pIdx=1:length(kinTrack.f)
        fIdx=kinTrack.f(pIdx);
        MTTipsAtTime=(EB3StartFrames<=fIdx)&(EB3EndFrames>=fIdx); % EB3 and Kin Co-exist if EB3 appear when Kin is alive.
        if(any(MTTipsAtTime))
            % associate EB3 to kin
            % Compute angle
            % threshold manually
            
            P1K=[kinTrack.x(pIdx)-P1P2Ref.origin(fIdx,1); ...
                kinTrack.y(pIdx)-P1P2Ref.origin(fIdx,2);  ...
                kinTrack.z(pIdx)-P1P2Ref.origin(fIdx,3)];
            
            P2K=[kinTrack.x(pIdx)-P2P1Ref.origin(fIdx,1); ...
                kinTrack.y(pIdx)-P2P1Ref.origin(fIdx,2);  ...
                kinTrack.z(pIdx)-P2P1Ref.origin(fIdx,3)];
            
            kinRhoP1=kinTrack.pole1.rho(pIdx);
            kinRhoP2=kinTrack.pole2.rho(pIdx);
            
            % Associate EB3 to a potential PKin Ref
            % Retrieve the associated Rho for each Pole
            P1KinAssociatedMTIndex=[];
            P2KinAssociatedMTIndex=[];
            MTTipsAtTimeIdx=find(MTTipsAtTime);
            EB3RhoP1=[];
            EB3RhoP2=[];
            EB3IdxP1=[];
            EB3IdxP2=[];
            for j=1:length(MTTipsAtTimeIdx)
                ebIdx=MTTipsAtTimeIdx(j);
                ebPIdx=(EB3Tracks(ebIdx).f==fIdx);
                aEB3RhoP1=EB3Tracks(ebIdx).pole1.rho(ebPIdx);
                aEB3RhoP2=EB3Tracks(ebIdx).pole2.rho(ebPIdx);
                EB3ClosestPoleStart=EB3Tracks(ebIdx).poleId(ebPIdx);
                if(aEB3RhoP1<kinRhoP1)&&(EB3ClosestPoleStart==1)
                    P1KinAssociatedMTIndex=[P1KinAssociatedMTIndex ebIdx];
                    EB3RhoP1=[EB3RhoP1 aEB3RhoP1];
                    EB3IdxP1=[EB3IdxP1 find(ebPIdx)];
                end
                if(aEB3RhoP2<kinRhoP2)&&(EB3ClosestPoleStart==2)
                    P2KinAssociatedMTIndex=[P2KinAssociatedMTIndex ebIdx];
                    EB3RhoP2=[EB3RhoP2 aEB3RhoP2];
                    EB3IdxP2=[EB3IdxP2 find(ebPIdx)];
                end
            end
            
            %% Retrieve full vector to estimate angle.
            MTVectorP1KinRef=cell2mat(arrayfun(@(t) [ t.x(t.f==fIdx)-P1P2Ref.origin(fIdx,1);  ...
                t.y(t.f==fIdx)-P1P2Ref.origin(fIdx,2);  ...
                t.z(t.f==fIdx)-P1P2Ref.origin(fIdx,3)], ...
                EB3Tracks(P1KinAssociatedMTIndex),'unif',0)');
            
            MTVectorP2KinRef=cell2mat(arrayfun(@(t) [ t.x(t.f==fIdx)-P2P1Ref.origin(fIdx,1);  ...
                t.y(t.f==fIdx)-P2P1Ref.origin(fIdx,2);  ...
                t.z(t.f==fIdx)-P2P1Ref.origin(fIdx,3)], ...
                EB3Tracks(P2KinAssociatedMTIndex),'unif',0)');
            
            %% Apply angle and distance cutoff
            MTAnglesP1Kin= vectorAngleND(MTVectorP1KinRef,P1K);
            MTAnglesP2Kin= vectorAngleND(MTVectorP2KinRef,P2K);
            angleCutoff=pi/2;
            if(strcmp(p.distType,'normalDist'))
                distP1=sin(MTAnglesP1Kin).*EB3RhoP1;
                associatedMTP1SubIdx=(((EB3RhoP1>p.minDistFromPole)&(distP1<cutoff)&(abs(MTAnglesP1Kin)<angleCutoff)));%
                distP2=sin(MTAnglesP2Kin).*EB3RhoP2;
                associatedMTP2SubIdx=(((EB3RhoP2>p.minDistFromPole)&(distP2<cutoff)&(abs(MTAnglesP2Kin)<angleCutoff)));%
                if(any(associatedMTP1SubIdx))
                    kinTrack.associatedTipsP1=[kinTrack.associatedTipsP1; EB3Tracks(P1KinAssociatedMTIndex(associatedMTP1SubIdx))] ;
                    kinTrack.associatedTipsP1Idx=[kinTrack.associatedTipsP1Idx EB3IdxP1(associatedMTP1SubIdx)];                   
                    kinTrack.distTipsP1=[kinTrack.distTipsP1, distP1(associatedMTP1SubIdx)] ;
                end
                if(any(associatedMTP2SubIdx))
                    kinTrack.associatedTipsP2=[kinTrack.associatedTipsP2; EB3Tracks(P2KinAssociatedMTIndex(associatedMTP2SubIdx))] ;
                    kinTrack.associatedTipsP2Idx=[kinTrack.associatedTipsP2Idx EB3IdxP2(associatedMTP2SubIdx)];                                      
                    kinTrack.distTipsP2=[kinTrack.distTipsP2, distP2(associatedMTP2SubIdx)] ;
                end
            else
                associatedMTP1=EB3Tracks(P1KinAssociatedMTIndex(abs(MTAnglesP1Kin)<cutoff));
                associatedMTP2=EB3Tracks(P2KinAssociatedMTIndex(abs(MTAnglesP2Kin)<cutoff))    ;
                if(~isempty(associatedMTP1))
                    kinTrack.associatedTipsP1=[kinTrack.associatedTipsP1; associatedMTP1] ;
                end
                if(~isempty(associatedMTP2))
                    kinTrack.associatedTipsP2=[kinTrack.associatedTipsP2; associatedMTP2] ;
                end
            end
            
            
        end
    end
end
end

function angle=vectorAngleND(a,b)
if(isempty(a))
    angle=[];
else
    if (size(a,2)==1)
        a=repmat(a,1,size(b,2));
    end
    if (size(b,2)==1)
        b=repmat(b,1,size(a,2));
    end
    f=@(a,b)atan2(norm(cross(a,b)), dot(a,b));
    angle=arrayfun(@(i) f(a(:,i),b(:,i)),1:size(a,2));
end
end
