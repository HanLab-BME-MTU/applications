function [manifoldCell,subManifoldCell]=buildFiberManifold(poles,kinTracks,tracks,distCutoff,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('kinDistCutoff',[])
ip.parse(varargin{:});
p=ip.Results;

manifoldCell=cell(length(poles),length(kinTracks));
subManifoldCell=cell(length(poles),length(kinTracks));

for kinIdx=1:length(kinTracks)
  kinTrack=kinTracks(kinIdx);
  progressText(kinIdx/length(kinTracks),'Build manifolds');
  for poIdx=1:length(poles)
    P=poles(poIdx);
    manifold=[P.getOverlapping(kinTrack) kinTrack];
    subManifold=manifold;
    if(~isempty(p.kinDistCutoff) )
      manifVector= manifold(2).getAddCoord(manifold(1).getMultCoord(-1));
      manifNorm=(manifVector.x.^2 + manifVector.y.^2 + manifVector.z.^2).^0.5;
      subManifold=[manifold(2).getAddCoord(manifVector.getMultCoord(p.kinDistCutoff(1)./manifNorm)) manifold(2).getAddCoord(manifVector.getMultCoord(p.kinDistCutoff(2)./manifNorm))];
    end
    manifoldCell{poIdx,kinIdx}=manifold;
    subManifoldCell{poIdx,kinIdx}=subManifold;
  end
end
