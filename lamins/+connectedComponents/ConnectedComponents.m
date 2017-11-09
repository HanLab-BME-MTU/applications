classdef ConnectedComponents < handle
    properties
        cc;
    end
    properties ( Transient )
        lm
    end
    properties ( Dependent )
        NumObjects
        ImageSize
        PixelIdxList
        Connectivity
    end
    methods
        function obj = ConnectedComponents(cc)
            obj.cc = cc;
        end
        function o = get.NumObjects(obj)
            o = obj.cc.NumObjects;
        end
        function o = get.ImageSize(obj)
            o = obj.cc.ImageSize;
        end
        function o = get.PixelIdxList(obj)
            o = obj.cc.PixelIdxList;
        end
        function o = get.Connectivity(obj)
            o = obj.cc.Connectivity;
        end
        function o = get.lm(obj)
            if(isempty(obj.lm))
                obj.lm = labelmatrix(obj.cc);
            end
            o = obj.lm;
        end
        function set.cc(obj,cc)
            obj.cc = cc;
            obj.lm = [];
        end
        function L = labelmatrix(obj)
            L = labelmatrix(obj.cc);
            obj.lm = L;
        end
        function rp = regionprops(obj,varargin)
            rp = regionprops(obj.cc,varargin{:});
        end
        function obj = dilate(obj,se)
            obj.lm = [];
            obj.cc = connectedComponents.ccDilate(obj.cc,se)
        end
        function obj = erode(obj,se)
            obj.lm = [];
            obj.cc = connectedComponents.ccErode(obj.cc,se);
        end
        function o = plus(a,b)
            import connectedComponents.*;
            o = ConnectedComponents(ccPlus(a.cc,b.cc));
        end 
        function o = minus(a,b)
            import connectedComponents.*;
            o = ConnectedComponents(ccMinus(a.cc,b.cc));
        end
        function o = not(obj)
            o = ~obj.lm;
        end
    end
end
