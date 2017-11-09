classdef Detections <  handle  & matlab.mixin.Copyable & dynamicprops
    % Data encapsulator for detections
    properties %(SetAccess = protected)
       xCoord;    % 3xN or 2xN coordinate
       yCoord;    % Image intensity
       zCoord;
       amp;       % Image amplitude (when estimated) 
       ref;       % optional FrameOfRef objects. 
    end
    
    methods
        
        function obj = Detections(movieInfo)
            if nargin>0
                for i=1:length(movieInfo)
                obj(i).xCoord=movieInfo(i).xCoord;
                obj(i).yCoord=movieInfo(i).yCoord;
                obj(i).zCoord=movieInfo(i).zCoord;
                obj(i).amp=movieInfo(i).amp;
                end
            end 
        end
        
%         function setFromPStruct(obj,aPstruct)
%             if isfield(aPstruct,'z')
%                 obj.xyz=[aPstruct.x; aPstruct.y; aPstruct.z]; 
%                 obj.dxyz=[aPstruct.x_pstd ; aPstruct.y_pstd; aPstruct.z_pstd];
%                 obj.scale=[aPstruct.s; aPstruct.s; aPstruct.s];
%             else 
%                 obj.xyz=[aPstruct.x; aPstruct.y;]; 
%                 obj.dxyz=[aPstruct.x_pstd ; aPstruct.y_pstd;];
%                 obj.scale=[aPstruct.s; aPstruct.s;];
%             end
%             obj.amp=aPstruct.A; obj.dAmp=aPstruct.A_pstd;
%             obj.bg=aPstruct.c; obj.dBg=aPstruct.c_pstd;
%         end
        
        function ret=is3D(obj)
            ret=~isempty(obj.zCoord);
        end     
        
        function pos=getPosMatrix(obj)
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.')
                det=obj(fIdx);

                pos=[det.xCoord(:,1),det.yCoord(:,1),det.zCoord(:,1)];
            end    
        end

        function pos=setFromSingleTrack(obj,track)
            for fIdx=1:length(obj)
                det=obj(fIdx);
                for i=1:length(track.x)
                    obj(track.f(i)).xCoord=[track.x(i)' 0.5*ones(length(track.lifetime,1))];
                    obj(track.f(i)).yCoord=[track.y(i)' 0.5*ones(length(track.lifetime,1))];
                    obj(track.f(i)).zCoord=[track.z(i)' 0.5*ones(length(track.lifetime,1))];
                end
            end    
        end

        function movieInfo=getStruct(obj)
            movieInfo(length(obj))=struct('xCoord',[],'yCoord',[],'zCoord',[]);
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.')
                det=obj(fIdx);
                movieInfo(fIdx).xCoord=det.xCoord;
                movieInfo(fIdx).yCoord=det.yCoord;
                movieInfo(fIdx).zCoord=det.zCoord;
                movieInfo(fIdx).amp=det.amp;
            end    
        end
          
          
        
        function addSphericalCoord(obj)
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                
                det=obj(fIdx);
                try
                    det.addprop('azimuth');      
                    det.addprop('elevation');    
                    det.addprop('rho');          
                catch
                end
                [det.azimuth,det.elevation,det.rho]=cart2sph(det.xCoord(:,1), ... 
                                                    det.yCoord(:,1), ...
                                                    det.zCoord(:,1));
    
            end    
        end
        
        function scale(obj,MD,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched=true;
            ip.addRequired('MD');
            ip.addOptional('isotropic',true);
            ip.parse(MD,varargin{:});
            p=ip.Results;
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                det=obj(fIdx);
                if(p.isotropic)
                    det.xCoord(:,1)=(det.xCoord(:,1)-1)*MD.pixelSize_+1;
                    det.yCoord(:,1)=(det.yCoord(:,1)-1)*MD.pixelSize_+1;
                    det.zCoord(:,1)=(det.zCoord(:,1)-1)*MD.pixelSize_+1;
                else
                    det.xCoord(:,1)=(det.xCoord(:,1)-1)*MD.pixelSize_+1;
                    det.yCoord(:,1)=(det.yCoord(:,1)-1)*MD.pixelSize_+1;
                    det.zCoord(:,1)=(det.zCoord(:,1)-1)*MD.pixelSizeZ_+1;
                end
                
            end
        end

        function addOffset(obj,X,Y,Z)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched=true;  
            p=ip.Results;
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                det=obj(fIdx);
                det.xCoord(:,1)=det.xCoord(:,1)+X;
                det.yCoord(:,1)=det.yCoord(:,1)+Y;
                det.zCoord(:,1)=det.zCoord(:,1)+Z;
            end
        end
        
        
        function obj=concatenate(obj,cdets)
            for fIdx=1:length(obj)
                det=obj(fIdx);
                det.xCoord=[det(fIdx).xCoord;cdets(fIdx).xCoord];
                det.yCoord=[det(fIdx).yCoord;cdets(fIdx).yCoord];
                det.zCoord=[det(fIdx).zCoord;cdets(fIdx).zCoord];
                det.amp=[det(fIdx).amp;cdets(fIdx).amp];
            end
        end
        
        function det=getSelectIdx(obj,indices)
            det=obj.copy();
            det.selectIdx(indices);
        end
        
        function obj=selectIdx(obj,indices)
            if(~iscell(indices))
                for fIdx=1:length(obj)
                    %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                    det=obj(fIdx);
                    det.xCoord=det.xCoord(indices,:);
                    det.yCoord=det.yCoord(indices,:);
                    det.zCoord=det.zCoord(indices,:);
                    det.amp=det.amp(indices,:);
                end
            else
                arrayfun(@(d,i) selectIdx(d,indices{i}),obj,1:length(indices),'unif',0);
            end
        end        
%         function setFromMovieInfo(obj,aMovieInfo)
%             if isfield(aMovieInfo,'zCoord')
%                 obj.xyz=[aMovieInfo.xCoord(:,1); aMovieInfo.yCoord(:,1); aMovieInfo.zCoord(:,1)]; 
%                 obj.dxyz=[aMovieInfo.xCoord(:,2); aMovieInfo.yCoord(:,2); aMovieInfo.zCoord(:,2)]; 
%             else 
%                 obj.xyz=[aMovieInfo.xCoord(:,1); aMovieInfo.yCoord(:,1)]; 
%                 obj.dxyz=[aMovieInfo.xCoord(:,2); aMovieInfo.yCoord(:,2)]; 
%             end
%             obj.amp=aMovieInfo.amp(:,1); obj.dAmp=aMovieInfo.amp(:,2);
%         end
%        
        
        function compare(obj)
        end 
        
        function setFromLabel(obj)
        end

%       function P = addprop(obj,propName)
%             if(~isscalar(obj))
%                 P = arrayfun(@(x) addprop(x,propName),obj,'UniformOutput',false);
%                 P = [P{:}];
%                 P = reshape(P,size(obj));
%             else
%                 P = addprop@dynamicprops(obj,propName);
%             end
%         end
%         
%         function out = subsasgn(A,S,B,varargin)
%             try
%                 if(isempty(A))
%                     A = Detections.empty;
%                 end
%                 out = builtin('subsasgn',A,S,B,varargin{:});
%             catch err
%                 switch(err.identifier)
%                     case 'MATLAB:noPublicFieldForClass'
%                         if(~all(isprop(A,S(1).subs)))
%                             rethrow(err);
%                         end
%                         if(nargin < 4)
%                             % Allow for [tracks.prop] = 5; for dynamic
%                             % properties
%                             out = arrayfun(@(t) subsasgn(t,S,B),A,'UniformOutput',false);
%                             out = [out{:}];
%                             out = reshape(out,size(A));
%                         else
%                             % Allow for
%                             % test = {1,2,3,4,5}
%                             % [tracks.prop] = test{:}
%                             % for dynamic properties
%                             out = arrayfun(@(t,b) subsasgn(t,S,b{1}),A,[{B} varargin],'UniformOutput',false);
%                             out = [out{:}];
%                             out = reshape(out,size(A));
%                         end
%                     otherwise
%                         rethrow(err)
%                 end
%             end
%         end
%         function varargout = subsref(A,S)
%             try
%                 [varargout{1:nargout}] = builtin('subsref',A,S);
%             catch err
%                 switch(err.identifier)
%                     case 'MATLAB:noSuchMethodOrField'
%                         if(all(isprop(A,S(1).subs)))
%                             % Allow for tracks.prop where prop is a dynamic
%                             % property
%                             varargout = arrayfun(@(t) subsref(t,S),A,'Unif',false);
%                         else
%                             rethrow(err);
%                         end
%                     otherwise
%                         rethrow(err);
%                 end
%             end
%         end
    end
end

