function mapDetectionsToKin(kinTracks,detections,cutoff,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('distType','normalDist');
ip.addParameter('minDistFromPole',500);
ip.addParameter('project',false);
ip.parse(varargin{:});
p=ip.Results;
 
% Pole info to circumvent projecion bug
P1P2Ref=kinTracks(1).pole1.ref;
P2P1Ref=kinTracks(1).pole2.ref;

for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'mapDetectionToKin.');
    try
        kinTrack.addprop('associatedDetectP1');
        kinTrack.addprop('associatedDetectP2');
    catch
    end;
    %        kinTrack.appearingMT=cell(1,length(kinTrack.poleRef));
    tmp(length(detections))=Detections();
    kinTrack.associatedDetectP1=tmp;
    tmp(length(detections))=Detections();
    kinTrack.associatedDetectP2=tmp;
    for pIdx=1:length(kinTrack.f)
        fIdx=kinTrack.f(pIdx);
        detect=detections(fIdx);
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
        P1KinAssociatedMTIndex=find(detect.pole1.rho<kinRhoP1);
        P2KinAssociatedMTIndex=find(detect.pole2.rho<kinRhoP2);
        
        EB3RhoP1=detect.pole1.rho(P1KinAssociatedMTIndex)';
        EB3RhoP2=detect.pole2.rho(P2KinAssociatedMTIndex)';

        
        %% Retrieve full vector to estimate angle.
        MTVectorP1KinRef=[  detect.xCoord(P1KinAssociatedMTIndex,1)-P1P2Ref.origin(fIdx,1), ...
                            detect.yCoord(P1KinAssociatedMTIndex,1)-P1P2Ref.origin(fIdx,2), ...
                            detect.zCoord(P1KinAssociatedMTIndex,1)-P1P2Ref.origin(fIdx,3) ]';
        
        MTVectorP2KinRef=[  detect.xCoord(P2KinAssociatedMTIndex,1)-P2P1Ref.origin(fIdx,1), ...
                            detect.yCoord(P2KinAssociatedMTIndex,1)-P2P1Ref.origin(fIdx,2), ...
                            detect.zCoord(P2KinAssociatedMTIndex,1)-P2P1Ref.origin(fIdx,3) ]';
        
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
                associatedDetectP1=detect.copy();
                associatedDetectP1.selectIdx(P1KinAssociatedMTIndex(associatedMTP1SubIdx));
                kinTrack.associatedDetectP1(fIdx)=associatedDetectP1 ;
            end
            if(any(associatedMTP2SubIdx))
                associatedDetectP2=detect.copy();
                associatedDetectP2.selectIdx(P2KinAssociatedMTIndex(associatedMTP2SubIdx));
                kinTrack.associatedDetectP2(fIdx)=associatedDetectP2 ;
            end
        else
%             associatedMTP1=detections(P1KinAssociatedMTIndex(abs(MTAnglesP1Kin)<cutoff));
%             associatedMTP2=detections(P2KinAssociatedMTIndex(abs(MTAnglesP2Kin)<cutoff))    ;
%             if(~isempty(associatedMTP1))
%                 kinTrack.associatedTipsP1=[kinTrack.associatedTipsP1; associatedMTP1] ;
%             end
%             if(~isempty(associatedMTP2))
%                 kinTrack.associatedTipsP2=[kinTrack.associatedTipsP2; associatedMTP2] ;
%             end
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
