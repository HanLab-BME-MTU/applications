function mapMTToKin(kinTracks,EB3Tracks,cutoff,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('distType','angle');
ip.addParameter('position','start')
ip.addParameter('kinDistCutoff',[])
ip.parse(varargin{:});
p=ip.Results;

EB3EndTime=[];
if(strcmp(p.position,'start'))
    EB3RhoP1=arrayfun(@(t) t.pole1.rho(1),EB3Tracks);
    EB3RhoP2=arrayfun(@(t) t.pole2.rho(1),EB3Tracks);
    EB3EndTime=[EB3Tracks.startFrame]';
else
    EB3RhoP1=arrayfun(@(t) t.pole1.rho(end),EB3Tracks);
    EB3RhoP2=arrayfun(@(t) t.pole2.rho(end),EB3Tracks);
    EB3EndTime=[EB3Tracks.endFrame]';
end;

EB3ClosestPoleStart=arrayfun(@(t) t.poleId(1),EB3Tracks);

for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'mapMTDisappearanceToKin.');
    try
        kinTrack.addprop('associatedMTP1');
        kinTrack.addprop('associatedMTP2');
    catch
    end;
    %        kinTrack.disappMT=cell(1,length(kinTrack.poleRef));
    kinTrack.associatedMTP1=[];
    kinTrack.associatedMTP2=[];
    for pIdx=1:length(kinTrack.f)
        fIdx=kinTrack.f(pIdx);
        coexistingEB3=(EB3EndTime==fIdx); % EB3 and Kin Co-exist if EB3 appear when Kin is alive.
        if(any(coexistingEB3))
            % associate EB3 to kin
            % Compute angle
            % threshold manually
            
            % Kinetochore distance from each pole at frame <fIdx>
            kinRhoP1=kinTrack.pole1.rho(pIdx);
            kinRhoP2=kinTrack.pole2.rho(pIdx);
            
            % Associate MT appearance to kinPole refererial using
            % Distance to pole
            P1KinAssociatedMT=(EB3RhoP1<kinRhoP1)&coexistingEB3&(EB3ClosestPoleStart==1);
            P2KinAssociatedMT=(EB3RhoP2<kinRhoP2)&coexistingEB3&(EB3ClosestPoleStart==2);
            if(~isempty(p.kinDistCutoff))
                P1KinAssociatedMT=P1KinAssociatedMT&(EB3RhoP1<(p.kinDistCutoff(2)+kinRhoP1))&(EB3RhoP1>(p.kinDistCutoff(1)+kinRhoP1));
                P2KinAssociatedMT=P2KinAssociatedMT&(EB3RhoP2<(p.kinDistCutoff(2)+kinRhoP2))&(EB3RhoP2>(p.kinDistCutoff(1)+kinRhoP2));
            end            
            P1KinAssociatedMTIndex=find(P1KinAssociatedMT);
            P2KinAssociatedMTIndex=find(P2KinAssociatedMT);
            
            %% Angle to kinPole axis
            if(strcmp(p.position,'start'))
                MTVectorP1KinRef=cell2mat(arrayfun(@(t) [t.pole1.x(1);t.pole1.y(1);t.pole1.z(1)],EB3Tracks(P1KinAssociatedMT),'unif',0)');
                MTVectorP2KinRef=cell2mat(arrayfun(@(t) [t.pole2.x(1);t.pole2.y(1);t.pole2.z(1)],EB3Tracks(P2KinAssociatedMT),'unif',0)');
            else
                MTVectorP1KinRef=cell2mat(arrayfun(@(t) [t.pole1.x(end);t.pole1.y(end);t.pole1.z(end)],EB3Tracks(P1KinAssociatedMT),'unif',0)');
                MTVectorP2KinRef=cell2mat(arrayfun(@(t) [t.pole2.x(end);t.pole2.y(end);t.pole2.z(end)],EB3Tracks(P2KinAssociatedMT),'unif',0)');
            end
            
            kinVectorP1KinRef=[kinTrack.pole1.x(pIdx); kinTrack.pole1.y(pIdx); kinTrack.pole1.z(pIdx) ];
            kinVectorP2KinRef=[kinTrack.pole2.x(pIdx); kinTrack.pole2.y(pIdx); kinTrack.pole2.z(pIdx) ];
            
            %%
            MTAnglesP1Kin= vectorAngleND(MTVectorP1KinRef,kinVectorP1KinRef);
            MTAnglesP2Kin= vectorAngleND(MTVectorP2KinRef,kinVectorP2KinRef);
            angleCutoff=0.2;
            switch p.distType
                case 'normalDistAndAngle'
                    distP1=sin(MTAnglesP1Kin).*EB3RhoP1(P1KinAssociatedMT)';
                    assocMTP1Kin=EB3Tracks(P1KinAssociatedMTIndex((distP1<cutoff)&(abs(MTAnglesP1Kin)<angleCutoff)));
                    distP2=sin(MTAnglesP2Kin).*EB3RhoP2(P2KinAssociatedMT)';
                    assocMTP2Kin=EB3Tracks(P2KinAssociatedMTIndex((distP2<cutoff)&(abs(MTAnglesP2Kin)<angleCutoff)));
                case 'normalDist'
                    distP1=sin(MTAnglesP1Kin).*EB3RhoP1(P1KinAssociatedMT)';
                    assocMTP1Kin=EB3Tracks(P1KinAssociatedMTIndex((distP1<cutoff)));
                    distP2=sin(MTAnglesP2Kin).*EB3RhoP2(P2KinAssociatedMT)';
                    assocMTP2Kin=EB3Tracks(P2KinAssociatedMTIndex((distP2<cutoff)));
                otherwise
                    assocMTP1Kin=EB3Tracks(P1KinAssociatedMTIndex(abs(MTAnglesP1Kin)<cutoff));
                    assocMTP2Kin=EB3Tracks(P2KinAssociatedMTIndex(abs(MTAnglesP2Kin)<cutoff));
            end
            
            kinTrack.associatedMTP1=[kinTrack.associatedMTP1; assocMTP1Kin] ;
            kinTrack.associatedMTP2=[kinTrack.associatedMTP2; assocMTP2Kin] ;
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
