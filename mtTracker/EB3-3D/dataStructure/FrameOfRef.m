classdef FrameOfRef < handle  & matlab.mixin.Copyable
   properties
      origin;
      X; % N*3 double
      Y; % N*3 double
      Z; % N*3 double
      frame; % N*1 uint16
   end
   methods
      function obj=setOriginFromTrack(obj,tr)
         obj.origin=[tr.x' tr.y' tr.z'];
         obj.frame=tr.f;

      end
      function obj=setZFromTrack(obj,tr)
         % Only the overlapping frame can be kept
         [F,idxTr,idxObj] = intersect(tr.f,obj.frame);
         obj.frame=F;
         obj.origin=obj.origin(idxObj,:);
         obj.Z=[tr.x(idxTr)' tr.y(idxTr)' tr.z(idxTr)']-obj.origin;
         obj.Z=obj.Z./repmat(sum(obj.Z.^2,2).^0.5,1,3);
      end

      function obj=genBaseFromZ(obj,trX)
          if(nargin<2)
            obj.X=[0*obj.Z(:,1),obj.Z(:,3),-obj.Z(:,2)];
            obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
          else
              [F,idxTr,idxObj] = intersect(trX.f,obj.frame);
              obj.frame=F;
              obj.X=[trX.x(idxTr)' trX.y(idxTr)' trX.z(idxTr)']-obj.origin;
              obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
              obj.X=obj.X-repmat(sum(obj.X.*obj.Z,2),1,3).*obj.Z;
              obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
          end
          obj.Y=cross(obj.X,obj.Z);
      end


      
      function obj=genCanonicalBase(obj)
          obj.X=repmat([1 0 0],[size(obj.origin,1) 1]);
          obj.Y=repmat([0 1 0],[size(obj.origin,1) 1]);
          obj.Z=repmat([0 0 1],[size(obj.origin,1) 1]);
      end
      function obj=genCanonicalRef(obj,frameNb)
          obj.origin=zeros(frameNb,3);
          obj.genCanonicalBase();
      end
      function tracks=getTracksFromBaseVector(obj,direction)
        vector=[]
        switch direction
        case 'X'
          vector=obj.X;
        case 'Y'
          vector=obj.Y;
        case 'Z'
          vector=obj.Z;
        otherwise
          error('Direction must be X,Y or Z')
        end
        tracks=TracksHandle();
        tracks.x=vector(:,1)';
        tracks.y=vector(:,2)';
        tracks.z=vector(:,3)';
        tracks.startFrame=obj.frame(1);
        tracks.endFrame=obj.frame(end);
      end
      function newBaseObject=applyBase(obj,tracksOrDetections,name)
        if(isa(tracksOrDetections,'Tracks'))
            newBaseObject=applyBaseToTrack(obj,tracksOrDetections,name);
        else
            newBaseObject=applyBaseToDetection(obj,tracksOrDetections,name);
        end
      end

      function tracksBase=applyBaseToTrack(obj,tracks,name)
          tracksBase=tracks.copy();
          try
              if(~(isempty(name)))
                  tracks.addprop(name);
              end
          catch
          end;
          for trIdx=1:length(tracks)
              tr=tracks(trIdx);
              trBase=tracksBase(trIdx);
              % Copying EB3 track
              trBase.addprop('ref');
              trBase.ref=obj;

              % Register in original tr
              try
                  trBase.addprop('originalRef');
              catch
              end;
              if(~(isempty(name)))
                  setfield(tr,name,trBase);
              end;
              trBase.originalRef=tr;
              for pIdx=1:length(tr.f)
                  f=tr.f(pIdx);
                  if(isempty(find(f==obj.frame,1)))
                      M=max(obj.frame);
                      if (f > M); f=M; else f=min(obj.frame); end;
                  end
                  B=obj.getBase(f);
                  orig=obj.getOrigAtFrame(f);
                  recentered=[(tr.x(pIdx)-orig(1)) (tr.y(pIdx)-orig(2)) (tr.z(pIdx)-orig(3))];
                  v=recentered*B;
                  trBase.x(pIdx)=v(1); trBase.y(pIdx)=v(2); trBase.z(pIdx)= v(3);
              end;
          end
        end
        function B= getBase(obj,f)
          pIdx=find(obj.frame==f,1);
          if(isempty(pIdx))
              if(f>max(obj.frame))   pIdx=length(obj.frame);  else   pIdx=1; end;
          end
          B=[obj.X(pIdx,:)' obj.Y(pIdx,:)' obj.Z(pIdx,:)'];
        end
      function orig= getOrigAtFrame(obj,f)
          pIdx=find(obj.frame==f,1);
          if(isempty(pIdx)); if(f<obj.frame); pIdx=1; else pIdx=length(obj.frame); end; end;
          orig=obj.origin(pIdx,:);
      end

      function detectionsBase=applyBaseToDetection(obj,detections,name)
          if(isempty(detections))
              detectionBase=[];
              return;
          end
          detectionsBase=detections.copy();
          try
              if(~(isempty(name)))
                  detections.addprop(name);
              end
          catch
          end;
          for fIdx=1:length(detections)
              detect=detections(fIdx);
              detectBase=detectionsBase(fIdx);
              % Copying EB3 track
              detectBase.ref=obj;
              f=min(fIdx,length(obj.X));
              B=obj.getBase(f);
              % easily optimized as implementend in poleDist
              for pIdx=1:size(detectBase.xCoord,1)
                  recentered=[(detect.xCoord(pIdx,1)-obj.origin(f,1)) (detect.yCoord(pIdx,1)-obj.origin(f,2)) (detect.zCoord(pIdx,1)-obj.origin(f,3))];
                  v=recentered*B;
                  detectBase.xCoord(pIdx,1)=v(1); detectBase.yCoord(pIdx,1)=v(2); detectBase.zCoord(pIdx,1)= v(3);
              end;
              if(~(isempty(name)))
                  setfield(detect,name,detectBase);
              end;
          end
      end
    end
end
