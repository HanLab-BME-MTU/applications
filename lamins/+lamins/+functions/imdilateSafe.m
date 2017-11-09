function I = imdilateSafe(I,se)
% imdilateSafe performs morphological dilation but only allows dilation into zero valued pixels.
    D = imdilate(I,se);
    I(~I) = D(~I);
end
