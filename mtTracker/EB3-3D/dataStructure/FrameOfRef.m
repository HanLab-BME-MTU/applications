classdef FrameOfRef < handle  & matlab.mixin.Copyable
   properties
      origin;
      X;
      Y;
      Z;
   end
   methods
      function setOriginFromTrack(obj,tr)
         obj.origin=[tr.x' tr.y' tr.z'];
      end
      function setZFromTrack(obj,tr)
         obj.Z=[tr.x' tr.y' tr.z']-obj.origin;
         obj.Z=obj.Z./repmat(sum(obj.Z.^2,2).^0.5,1,3);
      end
      
      function genBaseFromZ(obj)
          obj.X=[0*obj.Z(:,1),obj.Z(:,3),-obj.Z(:,2)];
          obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
          obj.Y=cross(obj.X,obj.Z);
      end
      
      function tracksBase=applyBaseToTrack(obj,tracks)
          tracksBase=tracks.copy();
          try
              tracks.addprop('altRef');
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
                  trBase.addprop('altRef');
                  tr.addprop('altRef');
              catch
              end;
              tr.altRef=[tr.altRef trBase];
              trBase.altRef=[trBase.altRef tr];
              for pIdx=1:length(tr.f)
                  f=min(tr.f(pIdx),length(obj.X));
                  B=[obj.X(f,:)' obj.Y(f,:)' obj.Z(f,:)'];
                  recentered=[(tr.x(pIdx)-obj.origin(f,1)) (tr.y(pIdx)-obj.origin(f,2)) (tr.z(pIdx)-obj.origin(f,3))];
                  v=recentered*B;
                  trBase.x(pIdx)=v(1); trBase.y(pIdx)=v(2); trBase.z(pIdx)= v(3);
              end;    
          end
      end
    end
end