classdef OrientationScaleSpaceResponse < OrientationSpaceResponse
    %OrientationScaleSpaceResponse Response object for
    %OrientationSpaceFilter at multiple scales
    
    properties
        scale
    end
    
    methods
        function obj = OrientationScaleSpaceResponse(filter,angularResponse,scale)
            obj = obj@OrientationSpaceResponse(filter,angularResponse);
            obj.scale = scale;
        end
        function Response = getResponseAtOrderFT(obj,K_new)
            Response = getResponseAtOrderFT@OrientationSpaceResponse(obj,K_new);
            Response = OrientationScaleSpaceResponse(Response.filter,Response.angularResponse,obj.scale);
        end
    end
    
end

