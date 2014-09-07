classdef ConnectedComponents < handle
    properties
        cc
        rp
    end
    methods
        function obj = ConnectedComponents(varargin)
            if(nargin > 0)
                if(islogical(varargin{1}))
                    obj.cc = bwconncomp(varargin{:});
                elseif(isstruct(varargin{1}))
                    obj.cc = varargin{1};
                    if(nargin > 1 && isstruct(varargin{2}))
                        obj.rp = varargin{2};
                    end
                end
            end
        end
        function filtered = filter(obj,criteria)
            newcc = filtercc(obj.cc,criteria);
            newrp = rp(criteria); 
            filtered = ConnectedComponents(newcc,newrp);
        end
        function rp = regionProps(obj,varargin);
            obj.rp = regionprops(obj.cc,varargin{:});
            rp = obj.rp;
        end
        function lm = labelMatrix(obj)
            lm = labelmatrix(obj.cc);
        end
        function show(obj)
            showPropMatrix(obj.cc);
        end
    end
end
