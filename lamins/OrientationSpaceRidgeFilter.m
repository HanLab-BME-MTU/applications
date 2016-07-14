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
            R(numel(obj)) = OrientationSpaceResponse;
            for o=1:numel(obj)
                R(o) = OrientationSpaceResponse(obj(o),angularResponse(:,:,:,o));
            end
            R = reshape(R,size(obj));
        end
    end
    
end

