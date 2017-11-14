function [mappedPoint,dist]=mapPointsTo1DManifold(points,manifold,cutoff,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('distType','normalDistPseudOptimized');
ip.addParameter('angleCutoff',[])
ip.parse(varargin{:});
p=ip.Results;
manifOrig=manifold(:,1);
manifVector=manifold(:,2)-manifOrig;
trackVector=points-repmat(manifOrig,1,size(points,2));


if(strcmp(p.distType,'normalDistPseudOptimized'))
    normManifVector=norm(manifVector);
    normalizeManifVector=manifVector/normManifVector;
    normVecRep=repmat(normalizeManifVector,1,size(points,2));
end

mappedPoint=[];
dist=[];
switch p.distType
    case 'normalDistAndAngle'
        %   distP1=sin(MTAnglesP1Kin).*norm(trackVector);
        %   assocMTP1Kin=TracksTracks(P1KinAssociatedMTIndex((distP1<cutoff)&(abs(tracksAngle)<p.angleCutoff)));
        error('not implemented')
    case 'normalDist'
        tracksAngle= vectorAngleND(trackVector,manifVector);
        tracksDist=sum(trackVector.^2,1).^(0.5);
        dist=sin(tracksAngle).*tracksDist;
        mappedPoint=((dist<cutoff)&(abs(tracksAngle)<pi/2)&(tracksDist<sum(manifVector.^2,1).^(0.5)));
    case 'euclideanDist'
        tracksAngle= vectorAngleND(trackVector,manifVector);
        tracksDist=sum(trackVector.^2,1).^(0.5);
        dist=sin(tracksAngle).*tracksDist;
        trackVector2=points-repmat(manifold(:,2),1,size(points,2));
        tracksDist2=sum(trackVector2.^2,1).^(0.5);
        mappedPoint=(((dist<cutoff)&(abs(tracksAngle)<pi/2)&(tracksDist<sum(manifVector.^2,1).^(0.5))) ... 
                    |(tracksDist<cutoff)|(tracksDist2<cutoff));
    case 'normalDistPseudOptimized'
        projPara=dot(trackVector,normVecRep);
        inBound=(projPara>0)&(projPara<normManifVector);
        mappedPoint=inBound;
        if(any(inBound))
            inBoundIdx=find(inBound);
            dist=sum((trackVector(:,inBound)- normVecRep(:,inBound).*repmat(projPara(inBound),3,1)).^2,1).^0.5;
            mappedPoint(inBoundIdx(dist>cutoff))=0;
        end
    case 'vertexDistOtsu'
        pointDist=points-repmat(manifold(:,1),1,size(points,2));
        minDist=sum(pointDist.^2,1).^(0.5);
        for mIdx=2:size(manifold,2)
            pointDist=points-repmat(manifold(:,mIdx),1,size(points,2));
            tmpMinDist=sum(pointDist.^2,1).^(0.5);
            minDist=min(tmpMinDist,minDist);
        end
        l=graythresh(double(minDist)/max(minDist))*max(minDist);
        mappedPoint=minDist<l;
    otherwise
        error('not implemented')
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
