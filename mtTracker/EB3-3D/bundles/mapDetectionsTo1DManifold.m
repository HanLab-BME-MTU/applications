function [mappedDetections]=mapDetectionsTo1DManifold(manifold,movieInfo,cutoff,varargin)
  % In this initial Q&D version, we use the angle with the vector to save some computation compared to the full blow reference
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addParameter('distType','normalDist');
  ip.addParameter('distCutoff',[])
  ip.addParameter('angleCutoff',[])
  ip.parse(varargin{:});
  p=ip.Results;

  % kinTrack.disappMT=cell(1,length(kinTrack.poleRef));
  mappedDetections(length(movieInfo))=struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);
  for pIdx=1:length(manifold(1).f)
      fIdx=manifold(1).f(pIdx);
      manifoldAtT=[[manifold(1).x(pIdx);manifold(1).y(pIdx);manifold(1).z(pIdx)], ...
          [manifold(2).x(pIdx);manifold(2).y(pIdx);manifold(2).z(pIdx) ]];
      points=[movieInfo(fIdx).xCoord(:,1),movieInfo(fIdx).yCoord(:,1),movieInfo(fIdx).zCoord(:,1)]';
      [mappedPoint,~]=mapPointsTo1DManifold(points,manifoldAtT,cutoff,varargin{:});
      mappedDetections(fIdx).xCoord=movieInfo(fIdx).xCoord(mappedPoint,:);
      mappedDetections(fIdx).yCoord=movieInfo(fIdx).yCoord(mappedPoint,:);
      mappedDetections(fIdx).zCoord=movieInfo(fIdx).zCoord(mappedPoint,:);      
      mappedDetections(fIdx).amp=movieInfo(fIdx).amp(mappedPoint,:);
      mappedDetections(fIdx).int=movieInfo(fIdx).int(mappedPoint,:);
  end
end

