classdef Skeleton < hgsetget
    properties
        parent
        edges
        vertices
    end
    methods
        function obj = Skeleton(image)
            if(nargin > 0)
                obj.parent = bwmorph(image,'skel',Inf);
            end
        end
        function get.edges
        end
        function get.vertices
        end
    end
end
