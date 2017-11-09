function [mask, thresh] = generateMask(Q,thresh)
    if(nargin < 2 || isempty(thresh))
        thresh = thresholdRosin(Q);
    end
    S = strel('disk',10);
    mask = imclose(Q > thresh,S);
    mask = imfill(mask,'holes');
    mask = imopen(mask,S);
end

%mask = imopen(imfill(imclose(Q > 1010,strel('disk',10)),'holes'),strel('disk',10));
