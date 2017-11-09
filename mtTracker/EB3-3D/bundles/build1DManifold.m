function [manifoldCell,subManifoldCell]=build1DManifold(trackSet1,trackSet2,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('kinDistCutoff',[])
ip.parse(varargin{:});
p=ip.Results;

manifoldCell=cell(length(trackSet1),length(trackSet2));
subManifoldCell=cell(length(trackSet1),length(trackSet2));

for kinIdx=1:length(trackSet2)
  kinTrack=trackSet2(kinIdx);
  progressText(kinIdx/length(trackSet2),'Build manifolds');
  for poIdx=1:length(trackSet1)
    P=trackSet1(poIdx);
    manifold=[P.getOverlapping(kinTrack) kinTrack];
    subManifold=manifold;
    if(~isempty(p.kinDistCutoff) )
      disp('creating subManifold')
      manifVector= manifold(2).getAddCoord(manifold(1).getMultCoord(-1));
      manifNorm=(manifVector.x.^2 + manifVector.y.^2 + manifVector.z.^2).^0.5;
      subManifold=[manifold(2).getAddCoord(manifVector.getMultCoord(p.kinDistCutoff(1)./manifNorm)) manifold(2).getAddCoord(manifVector.getMultCoord(p.kinDistCutoff(2)./manifNorm))];
    end
    manifoldCell{poIdx,kinIdx}=manifold;
    subManifoldCell{poIdx,kinIdx}=subManifold;
  end
end
