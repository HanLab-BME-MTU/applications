classdef projectDynROIProcess < ComputeMIPProcess & NonSingularProcess
    methods
        function obj = projectDynROIProcess(owner, varargin)
            obj = obj@ComputeMIPProcess(owner, varargin{:});
            obj=obj@NonSingularProcess(varargin{:});
            obj.funName_=(@(p) projectDynROI(owner,'processSingleProj',p));
        end
    end

    methods (Static)
        function name = getName()
            name = 'DynROI Maximum Intensity Projection';
        end
    end
end
