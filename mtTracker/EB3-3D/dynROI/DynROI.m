classdef DynROI < hgsetget & matlab.mixin.Copyable & handle

    properties (SetAccess = public, GetAccess = public)
        defaultRef; % If the ROI if 1D or 2D, it will generate a default frame of reference. 
    end

    methods (Abstract)
        [minmaxXBorder, minmaxYBorder,minmaxZBorder]=getBoundingBox(obj)
    end
end
