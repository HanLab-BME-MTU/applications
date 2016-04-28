classdef OrientationSpaceEdgeFilter < OrientationSpaceFilter
    %OrientationSpaceRidgeFilter OrientationSpace filter that does edge
    %filtering only
    
    properties
    end
    
    methods
        function obj = OrientationSpaceEdgeFilter(f_c,b_f,K)
            obj@OrientationSpaceFilter(f_c,b_f,K);
        end
        function R = getResponse(obj,I)
            If = fft2(I);
            angularResponse = obj.applyEdgeFilter(If);
            R = OrientationSpaceResponse(obj,angularResponse);
        end
    end
    
end

