function [ extended, extended_raw ] = extendVectorUntilConnected( bw, pts, vector, maxPixels, connectivity)
%extendVectorUntilConnected Extends a vector from a point in a binary image
% until it is connected. Used for extending segments of a partial skeleton
%
% bw is the black and white image describing the skeleton
%
% pts a p x 2 matrix describing the starting point [x y] or [c r]
%   for extension and is   usually an endpoint of a segment from bwmorph
%   Alternatively a vector of linear indices for bw
%   Alternatively a logical vector
%
% vector is either a vector with length p containing the vector angles
%   or a p x 2 matrix containing the change in [x y] or [c r]
%
% maxPixels is the maximum pixels by which extend the vectors
%
% connectivity is either a scalar value 4 or 8 or a matrix describing 
%   the neighborhood connection filter


% maxPixels = 16;

if(size(vector,2) == 1)
    % vector describes the angle theta
    vector = [cos(theta) sin(theta)];
end

if(size(pts,2) == 2)
    c = pts(:,1);
    r = pts(:,2);
elseif(isvector(pts))
    if(islogical(pts))
        pts = find(pts);
    end    
    [r,c] = ind2sub(size(bw),pts);
else
    [r,c] = find(pts);
end

% apply the 1 norm so that max(u,v) == 1
n = max(abs(vector),[],2);
vector = vector ./ repmat(n,[1 2]);
u = vector(:,1);
v = vector(:,2);

pixelSequence = 1:maxPixels;

% consider using intline or bresenham (e.g. Bresenham's line algorithm)
uu = round(u*pixelSequence);
vv = round(v*pixelSequence);
uuc = uu+repmat(c,1,maxPixels);
vvr = vv+repmat(r,1,maxPixels);

assert(~isscalar(connectivity) || connectivity == 4 || connectivity == 8,'Connectivity must either be 4, 8 or a filter');

neighborsFilter = ones(3);
neighborsFilter(2,2) = 0;

extended = double(bw);
keepExtending = ~isnan(u);
% extend_idx = sub2ind([1024 1024],vvr,uuc);
% extended_diag = bwmorph(threshed.lm > 0,'diag',Inf);
% % keepExtending(keepExtending) = any(extended_diag(extend_idx(keepExtending,:)),2);
% extend_neighbors = zeros(size(keepExtending));



% neighborsFilter(1,1) = 0;
% neighborsFilter(3,1) = 0;
% neighborsFilter(3,3) = 0;
% neighborsFilter(1,3) = 0;

for i=1:maxPixels
    keepExtending(vvr(:,i) < 1) = false;
    keepExtending(uuc(:,i) < 1) = false;
    keepExtending(vvr(:,i) > size(bw,1)) = false;
    keepExtending(uuc(:,i) > size(bw,2)) = false;
    extend_idx = sub2ind(size(bw),vvr(keepExtending,i),uuc(keepExtending,i));
    
    % 2015 05 18
    % If we bump directly into something while extending, stop
    keepExtending(keepExtending) = ~bw(extend_idx);
    extend_idx = extend_idx(~bw(extend_idx));
    
    extended(extend_idx) = 1;
    neighbors = imfilter(extended,neighborsFilter);
    extend_neighbors = zeros(size(keepExtending));
    extend_neighbors(keepExtending) = neighbors(extend_idx);
    keepExtending = keepExtending & extend_neighbors < 3;
end

extended_raw = extended;

extended = bwmorph(extended,'skel',Inf);
extended = bwmorph(extended,'spur',Inf);

%figure;
%imshowpair(extended,bw);

%figure;
%regions = bwlabel(~extended,4);
%regions(regions == 1) = 0;
%showPropMatrix(regions);

end

