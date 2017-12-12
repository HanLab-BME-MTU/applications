function [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj,XLimit,YLimit,ZLimit,fIdx,tracksInMask,myColormap,colorIndx,varargin)
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addOptional('cumulative',false);
  ip.parse(varargin{:});
  p=ip.Results;

  minXBorder=XLimit(1);
  maxXBorder=XLimit(2);
  minYBorder=YLimit(1);
  maxYBorder=YLimit(2);
  minZBorder=ZLimit(1);
  maxZBorder=ZLimit(2);

  %% Print tracks on the projections

  if(isempty(colorIndx))
%    colorIndx=floor(linspace(1,255,length(tracksInMask)));
     colorIndx=ones(1,length(tracksInMask));
  end
  
  if(isempty(myColormap))
      myColormap=255*jet(length(unique(colorIndx)));
  end

  if(~isempty(tracksInMask))
      keepIdx=false(1,length(tracksInMask));
        for tIdx=1:length(tracksInMask)
          tr=tracksInMask(tIdx);
          coordIdx=(tr.f==fIdx);
          if(any(coordIdx))
            X=tr.x(coordIdx);Y=tr.y(coordIdx); Z=tr.z(coordIdx);
            keepIdx(tIdx)= (Z>=minZBorder)&&(Z<=maxZBorder)&&(X>=minXBorder)&&(X<=maxXBorder)&&(Y>=minYBorder)&&(Y<=maxYBorder);
          end
      end
      tracksInMask=tracksInMask(keepIdx);
      colorIndx=colorIndx(keepIdx);
  end

  
  if(~isempty(tracksInMask))
    tracksXY=trackBinaryOverlay(XYProj,[minXBorder maxXBorder],[minYBorder maxYBorder],tracksInMask,fIdx,colorIndx,myColormap,varargin{:});
  else
    tracksXY=XYProj;
  end

  if(~isempty(tracksInMask))
      numTrack=length(tracksInMask);
      transTracks(numTrack)=struct('x',[],'y',[],'f',[]);
      for tIdx=1:numTrack
          transTracks(tIdx).x=tracksInMask(tIdx).z;
          transTracks(tIdx).y=tracksInMask(tIdx).y;
          transTracks(tIdx).f=tracksInMask(tIdx).f;
      end

%     capturedEB3ZY=tracksInMask.copy();
%     for tIdx=1:length(capturedEB3ZY)
%       capturedEB3ZY(tIdx).x=tracksInMask(tIdx).z ;%*MD.pixelSize_/MD.pixelSizeZ_;
%     end
    tracksZY=trackBinaryOverlay(ZYProj,[minZBorder maxZBorder],[minYBorder maxYBorder],transTracks,fIdx,colorIndx,myColormap,varargin{:});
  else
    tracksZY=ZYProj;
  end;


  if(~isempty(tracksInMask))
      for tIdx=1:numTrack
          transTracks(tIdx).y=tracksInMask(tIdx).x;
      end

%     capturedEB3ZX=capturedEB3ZY.copy();
%     for tIdx=1:length(capturedEB3ZX)
%       capturedEB3ZX(tIdx).y=tracksInMask(tIdx).x;
%     end
    tracksZX=trackBinaryOverlay(ZXProj,[minZBorder maxZBorder],[minXBorder maxXBorder],transTracks,fIdx,colorIndx,myColormap,varargin{:});
  else
    tracksZX=ZXProj;
  end;
