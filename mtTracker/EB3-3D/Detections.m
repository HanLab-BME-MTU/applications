classdef Detections < hgsetget
    % Data encapsulator for detections
    properties %(SetAccess = protected)
       xyz; dxyz;    % 3xN or 2xN coordinate
       I; dI;        % Image intensity
       amp; dAmp;        % Image amplitude (when estimated) 
       scale; dScale;        % 3xN or 2xN scales
       bg; dBg;
       label;        % An image representing the label of the proposed object.
       MD;           % associated MovieData
       frameIdx;     % Idx in the associated movieData
    end
    
    methods
        
        function obj = Detections(aMD,aFrameIdx)
            if nargin>0
            obj.MD=aMD;
            obj.frameIdx=aFrameIdx;
            end 
        end
        
        function setFromPStruct(obj,aPstruct)
            if isfield(aPstruct,'z')
                obj.xyz=[aPstruct.x; aPstruct.y; aPstruct.z]; 
                obj.dxyz=[aPstruct.x_pstd ; aPstruct.y_pstd; aPstruct.z_pstd];
                obj.scale=[aPstruct.s; aPstruct.s; aPstruct.s];
            else 
                obj.xyz=[aPstruct.x; aPstruct.y;]; 
                obj.dxyz=[aPstruct.x_pstd ; aPstruct.y_pstd;];
                obj.scale=[aPstruct.s; aPstruct.s;];
            end
            obj.amp=aPstruct.A; obj.dAmp=aPstruct.A_pstd;
            obj.bg=aPstruct.c; obj.dBg=aPstruct.c_pstd;
        end
        
        function answer=is3D(obj)
            answer=size(obj.xyz,2)==3;
        end     
        
        function setFromMovieInfo(obj,aMovieInfo)
            if isfield(aMovieInfo,'zCoord')
                obj.xyz=[aMovieInfo.xCoord(:,1); aMovieInfo.yCoord(:,1); aMovieInfo.zCoord(:,1)]; 
                obj.dxyz=[aMovieInfo.xCoord(:,2); aMovieInfo.yCoord(:,2); aMovieInfo.zCoord(:,2)]; 
            else 
                obj.xyz=[aMovieInfo.xCoord(:,1); aMovieInfo.yCoord(:,1)]; 
                obj.dxyz=[aMovieInfo.xCoord(:,2); aMovieInfo.yCoord(:,2)]; 
            end
            obj.amp=aMovieInfo.amp(:,1); obj.dAmp=aMovieInfo.amp(:,2);
        end
        
        
        function getMovieInfo(obj)
            
        end
        
        function compare(obj)
        end 
        
        function setFromLabel(obj)
        end
                                    
        %function setFromLabel(label)
%         function showLabel(vol)
%             ip=inputParser;
%             ip.add
%     
%         end
    end
end

