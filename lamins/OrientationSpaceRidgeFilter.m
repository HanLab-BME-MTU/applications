classdef OrientationSpaceRidgeFilter < OrientationSpaceFilter
    %OrientationSpaceRidgeFilter OrientationSpace filter that does line
    %filtering only
    
    properties
    end
    
    methods
        function obj = OrientationSpaceRidgeFilter(f_c,b_f,K)
            obj@OrientationSpaceFilter(f_c,b_f,K);
        end
        function R = getResponse(obj,I)
            If = fft2(I);
            angularResponse = obj.applyRidgeFilter(If);
            R = OrientationSpaceResponse(obj,angularResponse);
        end
    end
    
end

