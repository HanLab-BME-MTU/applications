function [tracksXY,tracksZY,tracksZX]=overlayProjDetections(XYProj,ZYProj,ZXProj,XLimit,YLimit,ZLimit,detections,myColormap,colorIndx,varargin)
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
    colorIndx=arrayfun(@(d) ones(1,length(d.zCoord(:,1))),detections,'unif',0);
  end

  if(~iscell(colorIndx))
    colorIndx={colorIndx};
  end
  
  if(isempty(myColormap))
    myColormap=255*jet(length(unique(vertcat(colorIndx{:}))));
  end
  
  if(~isempty(detections))
      keepIndx=cell(1,length(detections));
      for dIdx=1:length(detections)
          d=detections(dIdx);
          if(~isempty(d.xCoord))
              keepIndx{dIdx}=(d.zCoord(:,1)>minZBorder)&(d.zCoord(:,1)<maxZBorder)& ...
                  (d.xCoord(:,1)>minXBorder)&(d.xCoord(:,1)<maxXBorder)& ...
                  (d.yCoord(:,1)>minYBorder)&(d.yCoord(:,1)<maxYBorder);
          end
          colorIndx{dIdx}=colorIndx{dIdx}(keepIndx{dIdx});
      end
      detections=detections.copy().selectIdx(keepIndx);
  end


  if(~isempty(detections))
      % Only Keep detections within ZLimit
    tracksXY=detectionBinaryOverlay(XYProj,[minXBorder maxXBorder],[minYBorder maxYBorder],detections,colorIndx,myColormap,varargin{:});
  else
    tracksXY=XYProj;
  end

  if(~isempty(detections))
      numDet=length(detections);
      trdetections(numDet)=struct('xCoord',[],'yCoord',[]);
      for tIdx=1:numDet
          trdetections(tIdx).xCoord=detections(tIdx).zCoord;
          trdetections(tIdx).yCoord=detections(tIdx).yCoord;
      end
    
%     capturedEB3ZY=tracksInMask.copy();
%     for tIdx=1:length(capturedEB3ZY)
%       capturedEB3ZY(tIdx).x=tracksInMask(tIdx).z ;%*MD.pixelSize_/MD.pixelSizeZ_;
%     end
    tracksZY=detectionBinaryOverlay(ZYProj,[minZBorder maxZBorder],[minYBorder maxYBorder],trdetections,colorIndx,myColormap,varargin{:});
  else
    tracksZY=ZYProj;
  end;


  if(~isempty(detections))
      for tIdx=1:numDet
          trdetections(tIdx).yCoord=detections(tIdx).xCoord;
      end

%     capturedEB3ZX=capturedEB3ZY.copy();
%     for tIdx=1:length(capturedEB3ZX)
%       capturedEB3ZX(tIdx).y=tracksInMask(tIdx).x;
%     end
    tracksZX=detectionBinaryOverlay(ZXProj,[minZBorder maxZBorder],[minXBorder maxXBorder],trdetections,colorIndx,myColormap,varargin{:});
  else
    tracksZX=ZXProj;
  end;
