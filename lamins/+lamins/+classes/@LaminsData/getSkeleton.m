function S = getSkeleton(obj,varargin)
    I = obj.cellReader{varargin{:}};
    S = getLaminSkeleton(I);
end
