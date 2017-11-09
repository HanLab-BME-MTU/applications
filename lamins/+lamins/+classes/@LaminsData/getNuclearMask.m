function [mask,thresh] = getNuclearMask(obj,varargin)
    I = horzcat(obj.cellReader{varargin{:}});
    [mask, thresh] = generateMask(I);
end
