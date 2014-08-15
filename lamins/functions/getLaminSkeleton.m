function [ S, others ] = getLaminSkeleton( I, sigma )
%getLaminSkeleton Skeletonizes a lamin skeleton
if(nargin < 2)
    sigma = 4.5;
end

[res, theta, nms] = steerableDetector(double(I),2,sigma);
thresh = thresholdOtsu(res);
mask = generateMask(I);
B = mask & res > thresh;
S = bwmorph(B,'thin',Inf);

if(nargout > 1)
    others.res = res;
    others.theta = theta;
    others.nms = nms;
    others.mask = mask;
    others.B = B;
end

end

