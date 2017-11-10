function [ onion, onionLabeled ] = onionLabel( cc ,  initialIdx, adjacency)
%onionLabel Label connected components by connected distance from a
%specified connected component

if(nargin < 2)
    initialIdx = 1;
end

if(nargin < 3)
    adjacency = zeros(cc.NumObjects,cc.NumObjects);
    lm = labelmatrix(cc);
    dilated = imdilate(lm,strel('disk',3));
    lm(~lm) = dilated(~lm);
    nz_sub2ind = @(x,y) sub2ind(size(adjacency),x(x > 0 & y > 0),y(x > 0 & y > 0));
    % four connectivity
    adjacency(nz_sub2ind(lm(1:end-1,:),lm(2:end,:))) = 1;
    adjacency(nz_sub2ind(lm(:,1:end-1),lm(:,2:end))) = 1;
    % eight connectivity
    adjacency(nz_sub2ind(lm(1:end-1,1:end-1),lm(2:end,2:end))) = 1;
    adjacency(nz_sub2ind(lm(1:end-1,2:end),lm(2:end,1:end-1))) = 1;
    adjacency = adjacency | adjacency';
    adjacency = adjacency & ~eye(size(adjacency));
end

onion = zeros(1,cc.NumObjects);

adjacency = double(adjacency);
A = adjacency;
l = 1;

while(~all(onion) && l <= cc.NumObjects)
    candidates = A(initialIdx,:);
    onion(any(candidates,1) & ~onion) = l;
    l = l + 1
    A = A*adjacency;
end

onion(initialIdx) = 0;

if(nargout > 1)
    onionLabeled = lamins.functions.propmatrix(cc,onion);
end


end

