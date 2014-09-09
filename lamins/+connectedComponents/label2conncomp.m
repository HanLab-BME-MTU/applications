function [ out ] = label2conncomp( lm, conn)
%labelconncomp Converts label matrix back to connected components map
% lm is the label matrix from labelmatrix
% conn is the connectivity , which is ignored. Assumes 8

if(nargin < 2)
    conn = 8;
end
out.Connectivity = conn;
out.ImageSize = size(lm);
out.NumObjects = double(max(lm(:)));
% out.PixelIdxList = cell(1,out.NumObjects);
% 
% for i= 1:out.ImageSize
%     out.PixelIdxList{i} = find(lm == i);
% end
nzList = find(lm);
out.PixelIdxList = accumarray(lm(nzList),nzList,[],@(x) {x});


end

