classdef FrameOfRef < handle  & matlab.mixin.Copyable
   properties
      origin;
      X;
      Y;
      Z;
   end
   methods
      function obj=setOriginFromTrack(obj,tr)
         obj.origin=[tr.x' tr.y' tr.z'];
      end
      function obj=setZFromTrack(obj,tr)
         obj.Z=[tr.x' tr.y' tr.z']-obj.origin;
         obj.Z=obj.Z./repmat(sum(obj.Z.^2,2).^0.5,1,3);
      end

      function obj=genBaseFromZ(obj)
          obj.X=[0*obj.Z(:,1),obj.Z(:,3),-obj.Z(:,2)];
          obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
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
                  f=min(tr.f(pIdx),length(obj.X));
                  B=obj.getBase(f);
                  recentered=[(tr.x(pIdx)-obj.origin(f,1)) (tr.y(pIdx)-obj.origin(f,2)) (tr.z(pIdx)-obj.origin(f,3))];
                  v=recentered*B;
                  trBase.x(pIdx)=v(1); trBase.y(pIdx)=v(2); trBase.z(pIdx)= v(3);
              end;
          end
      end
      function B= getBase(obj,f)
             B=[obj.X(f,:)' obj.Y(f,:)' obj.Z(f,:)'];
      end
      function detectionsBase=applyBaseToDetection(obj,detections,name)
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
