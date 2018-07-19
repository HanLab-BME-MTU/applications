function [ out ] = label2conncomp( lm, conn)
%labelconncomp Converts label matrix back to connected components map
% lm is the label matrix from labelmatrix
% conn is the connectivity , which is ignored. Assumes 8

% Mark Kittisopikul
% UT Southwestern
% 2014 / 11 / 17

if(nargin < 2)
    conn = 8;
end
out.Connectivity = conn;
out.ImageSize = size(lm);
out.NumObjects = double(max(lm(:)));
% find nonzero values
nzList = find(lm);
% group label indices together
out.PixelIdxList = accumarray(lm(nzList),nzList,[],@(x) {x})';



end

