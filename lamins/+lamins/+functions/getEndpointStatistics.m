function [ stats ] = getEndpointStatistics( idx, cc, fits, nOrientations)
%getEndpointStatistics Summary of this function goes here
%   Detailed explanation goes here

[r, c] = ind2sub(cc.ImageSize,idx);
lm = labelmatrix(cc);
label = lm(idx);
% copy matrix of same size as lm
orderMap = lm;
for i=1:cc.NumObjects
    orderMap(cc.PixelIdxList{i}) = 1:length(cc.PixelIdxList{i});
end
order = orderMap(idx);

numAngles = (cellfun(@length,fits)-2)/2;

stats = [];

for i=1:min(max(numAngles),nOrientations)
    valid = numAngles >= i;
    orientation = cellfun(@(x) x(2*i-1),fits(valid));
    response = cellfun(@(x) x(2*i),fits(valid));
    stats = [stats ; [r(valid) c(valid) double(label(valid)) double(order(valid)) orientation response] ];
end

end

