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
            R(numel(obj)) = OrientationSpaceResponse;
            for o=1:numel(obj)
                R(o) = OrientationSpaceResponse(obj(o),angularResponse(:,:,:,o));
            end
            R = reshape(R,size(obj));
        end
    end
    
end

