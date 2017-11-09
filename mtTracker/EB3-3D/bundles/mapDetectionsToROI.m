function mappedDetections=mapDetectionsToROI(detections,roi,distCutoff,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('distType','vertexDistOtsu');
ip.addParameter('mappedTracksField','mapped')
ip.parse(varargin{:});
p=ip.Results;

mappedDetections=cell(1,length(detections));
for fIdx=1:length(detections)
    points=[detections(fIdx).xCoord(:,1)' ; detections(fIdx).yCoord(:,1)'; detections(fIdx).zCoord(:,1)'];
    pIdx=find(roi(1).f==fIdx);
    manifoldAtT=[[roi(1).x(pIdx);roi(1).y(pIdx);roi(1).z(pIdx)], ...
        [roi(2).x(pIdx);roi(2).y(pIdx);roi(2).z(pIdx) ]];
    [mappedPoint,~]=mapPointsTo1DManifold(points,manifoldAtT,distCutoff,'distType',p.distType);
    mappedDetections{fIdx}=mappedPoint;
end;

