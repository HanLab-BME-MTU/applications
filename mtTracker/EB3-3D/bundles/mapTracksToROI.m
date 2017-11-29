function [manifoldCell,subManifoldCell]=mapTracksToROI(tracks,ROIs,distCutoff,varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('distType','angle');
ip.addParameter('position','end')
ip.addParameter('manifoldEntry',true)
ip.addParameter('kinDistCutoff',[])
ip.addParameter('mappedTracksField','associatedMT')
ip.parse(varargin{:});
p=ip.Results;

manifoldCell=cell(length(poles),length(kinTracks));
subManifoldCell=cell(length(poles),length(kinTracks));
mappedTrackCell=cell(1,length(kinTracks));
parfor kinIdx=1:length(kinTracks)
  kinTrack=kinTracks(kinIdx);
  allMappedTracks=[];
  for poIdx=1:length(poles)
    P=poles(poIdx);
    manifold=[P.getOverlapping(kinTrack) kinTrack];
    subManifold=manifold;
    if(~isempty(p.kinDistCutoff) )
      manifVector= manifold(2).getAddCoord(manifold(1).getMultCoord(-1));
      manifNorm=(manifVector.x.^2 + manifVector.y.^2 + manifVector.z.^2).^0.5;
      subManifold=[manifold(2).getAddCoord(manifVector.getMultCoord(p.kinDistCutoff(1)./manifNorm)) manifold(2).getAddCoord(manifVector.getMultCoord(p.kinDistCutoff(2)./manifNorm))];
    end
    mappedTracks=mapTracksTo1DManifold(subManifold,tracks,distCutoff,'position',p.position);
    %% only keep MT that are aligned with the tube.
    % for each mapped track, take the last non-mapped point,
    % measure if the last unmap solution Z and XY are correct.
    % elevation i

    if(p.manifoldEntry)
        ref=buildRefsFromTracks(P,kinTrack);
        mappedTracksRef=ref.applyBase(mappedTracks,'');
        subManifoldREf=ref.applyBase(subManifold,'');
        manifoldEntry=false(1,length(mappedTracksRef));
        for tIdx=1:length(mappedTracksRef)
            % remap each timepoint to manifold
            mt=mappedTracksRef(tIdx);
            [~,manifSubIndx]=intersect(subManifoldREf(1).f,mt.f);
            nonMappedPoint=((((mt.x.^2+mt.y.^2).^(0.5))>distCutoff)| ...
                            (mt.z<subManifoldREf(1).z(manifSubIndx))| ...
                            (mt.z>subManifoldREf(2).z(manifSubIndx)));
            lastNonMapped=find(nonMappedPoint,1,'last');
            if(isempty(lastNonMapped))
              manifoldEntry(tIdx)=true;
            else
              manifoldEntry(tIdx)=((mt.x(lastNonMapped).^2+mt.y(lastNonMapped).^2).^0.5)<distCutoff;
            end;
        end
        mappedTracks=mappedTracks(manifoldEntry);
    end
    allMappedTracks=[allMappedTracks mappedTracks];
%     manifoldCell{poIdx,kinIdx}=manifold;
%     subManifoldCell{poIdx,kinIdx}=subManifold;
  end
  mappedTrackCell{kinIdx}=allMappedTracks;
end
for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);

  try
    kinTrack.addprop(p.mappedTracksField);
  catch
  end;
  setfield(kinTrack,p.mappedTracksField,mappedTrackCell{kinIdx});
end
