function mappedTracks=mapTracksTo1DManifold(manifold,tracks,cutoff,varargin)
  % In this initial Q&D version, we use the angle with the vector to save some computation compared to the full blow reference
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addParameter('distType','normalDist');
  ip.addParameter('position','start')
  ip.addParameter('distCutoff',[])
  ip.addParameter('angleCutoff',[])
  ip.parse(varargin{:});
  p=ip.Results;

  if(strcmp(p.position,'start'))
    TracksEndTime=[tracks.startFrame]';
  else
    TracksEndTime=[tracks.endFrame]';
  end;

  %        kinTrack.disappMT=cell(1,length(kinTrack.poleRef));
  mappedTracks=[];
  for pIdx=1:length(manifold(1).f)
    fIdx=manifold(1).f(pIdx);
    % associate tracks to manifold through time
    coexistingTracks=(TracksEndTime==fIdx); % Tracks and Kin Co-exist if Tracks appear when Kin is alive.
    if(any(coexistingTracks))
      % collect points
      manifoldPreAssociatonIndx=find(coexistingTracks);
      trackPoint=zeros(3,length(manifoldPreAssociatonIndx));
      for tIdxIdx=1:length(manifoldPreAssociatonIndx)
        tr=tracks(manifoldPreAssociatonIndx(tIdxIdx));
        if(strcmp(p.position,'start'))
          trackPoint(:,tIdxIdx)=[tr.x(1);tr.y(1);tr.z(1)];
        else
          trackPoint(:,tIdxIdx)=[tr.x(end);tr.y(end);tr.z(end)];
        end
      end
      manifoldAtT=[[manifold(1).x(pIdx);manifold(1).y(pIdx);manifold(1).z(pIdx)], ...
                   [manifold(2).x(pIdx);manifold(2).y(pIdx);manifold(2).z(pIdx) ]];
      [mappedPoint,~]=mapPointsTo1DManifold(trackPoint,manifoldAtT,cutoff,varargin{:});
      mappedTracks=[mappedTracks; tracks(manifoldPreAssociatonIndx(mappedPoint))];
    end
  end
end

