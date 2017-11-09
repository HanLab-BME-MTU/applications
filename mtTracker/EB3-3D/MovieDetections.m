classdef MovieDetections < matlab.mixin.SetGet
    % Data encapsulator for detections in a movie
    
    properties (SetAccess = protected)
       detectionsList;      
    end
    
    methods
        %% Constructor
        function obj = MovieDetections(varargin)
            ip=inputParser;
            ip.addOptional('detectionsCell',{});
            ip.parse(varargin{:})
            obj.detectionsList=ip.Results.detectionsCell;
        end


        %% convertors
        function setFromPStruct(obj,aPstruct)
            obj.detectionsList=cell(1,length(aPstruct));
            for i=1:length(aPstruct)
                det=Detections();
                set(det,'frameIdx',i);
                det.setFromPStruct(aPstruct(i));
                obj.detectionsList{i}=det;
            end
        end
    end
    methods (Static,Abstract)
    end
end

