%% Script settings
arrow_scale = 5;

%% Load initial data
if(~exist('Lamins','var'))
    initLamins;
end
% I = Lamins(2).cellReader{3,1,15};
% I = Lamins(2).cellReader{2,1,15};
I = im2double(imadjust(I,stretchlim(I,0)));

%% Do steerable detection
[res, theta, nms] = steerableDetector(double(I),4,5);
mask = Lamins(2).getNuclearMask(3,1,15);
mask = ones(size(I));
skel = nms ~= 0 & mask;
skel = bwmorph(skel,'skel',Inf);

%% Do thresholding by area
theta_stats = getSegmentOrientationStats(theta,skel);
areaThresh = thresholdRosin([theta_stats.rp.Area]);
areaThresh = 5;
threshed.cc = filtercc(theta_stats.cc,[theta_stats.rp.Area] > areaThresh);
threshed.lm = labelmatrix(threshed.cc);

%% Identify end points and branch points
endpts = bwmorph(threshed.lm,'endpoints');
branchpts = bwmorph(threshed.lm,'branchpoints');

%% Determine angle of segment near end points
% The dilation method would fail with nearby points
% labeled_endpts = threshed.lm;
% labeled_endpts(~endpts) = 0;
% dilated_labeled_endpts = imdilate(labeled_endpts,strel('disk',5));
% endpt_neighborhood_mask = dilated_labeled_endpts == threshed.lm;
endpt_neighborhood_mask = bwdistgeodesic(threshed.lm > 0,endpts) < 6;
endpt_neighborhood = threshed.lm;
endpt_neighborhood(~endpt_neighborhood_mask) = 0;
endpt_neighborhood_theta = getSegmentOrientationStats(theta,endpt_neighborhood > 0);
endpt_theta_matrix = propmatrix(endpt_neighborhood_theta.cc,[endpt_neighborhood_theta.rp.MeanIntensity]);

rpcentroid = regionprops(endpt_neighborhood_theta.cc,'Centroid');
centroids = vertcat(rpcentroid.Centroid);
centroid_idx = sub2ind([1024 1024],round(centroids(:,1)),round(centroids(:,2)));
centroid_matrix = propmatrix(endpt_neighborhood_theta.cc,centroid_idx);
[centroid.c,centroid.r] = ind2sub([1024 1024],centroid_matrix(find(endpts)));
% find and centroid coordinates are reversed
[endpt.r,endpt.c] = find(endpts);

%% Show matrix
showLabelMatrix(threshed.lm);

%% Show end points and centroids
hold on;
scatter(endpt.c,endpt.r)
scatter(centroid.c,centroid.r);
hold off;

%% Show centroid to endpoint vectors
hold on;
local.u = endpt.c - centroid.c;
local.v = endpt.r - centroid.r;
local.n = sqrt(local.u.^2 +local.v.^2)./arrow_scale;
local.x = endpt.c;
local.y = endpt.r;
% red arrows
quiver(endpt.c,endpt.r, local.u ./ local.n, local.v./ local.n ,0,'red');
hold off;

%% Show steerable filter orientation at endpoints
hold on
endpt_theta_matrix(~endpts)  = NaN;
%white arrows
[h, filter.x, filter.y, filter.u, filter.v ] = drawOrientation(endpt_theta_matrix);
filter.n = sqrt(filter.u.^2 +filter.v.^2);
hold off;

%% Breakup skeleton into segments and determine vectors between vertices
branchpt.idx = find(branchpts);
endpt.idx = find(endpts);

vertices = branchpts | endpts;
[edges_cc, vertices_cc, pairs] = bwtrace(threshed.lm > 0,vertices);

vlm = labelmatrix(vertices_cc);
vfilt = cellfun(@length,vertices_cc.edges(vlm(endpt.idx))) == 1;
v_endpts = vlm(endpt.idx(vfilt));

edge_vertices = [edges_cc.vertices{[vertices_cc.edges{v_endpts}]'}];
% reorder the vertices to have end points last
edge_vertices = [edge_vertices(edge_vertices ~= repmat(v_endpts',2,1)) v_endpts];
edge_vertices = reshape([vertices_cc.PixelIdxList{edge_vertices(:)}],[],2);

[seg.row, seg.col ] = ind2sub(vertices_cc.ImageSize,edge_vertices);

%% Plot segments
hold on;
seg.u = seg.col(:,2)-seg.col(:,1);
seg.v = seg.row(:,2)-seg.row(:,1);
seg.n = sqrt(seg.u.^2 +seg.v.^2)./arrow_scale;
quiver(seg.col(:,2),seg.row(:,2),seg.u ./ seg.n,seg.v ./ seg.n,0,'magenta');
hold off;

%% Initial segment extension along centroid to end point
u = endpt.c - centroid.c;
v = endpt.r - centroid.r;
n = sqrt(u.^2 +v.^2);
un = u ./ n;
vn = v ./ n;
extension_distance = 16;
uu = round(un*(1:extension_distance));
vv = round(vn*(1:extension_distance));
uuc = uu+repmat(endpt.c,1,extension_distance);
vvr = vv+repmat(endpt.r,1,extension_distance);
extend_I = zeros(1024);
for i=1:extension_distance
    extend_idx = sub2ind([1024 1024],vvr(:,i),uuc(:,i));
    extend_idx = extend_idx(~isnan(extend_idx));
    extend_I(extend_idx) = extension_distance - i + 1;
end

%% Show extension
figure;
imshowpair(extend_I,threshed.lm > 0);

%% Get labeled endpts
extension_distance = 16;
uu = round(un*(1:extension_distance));
vv = round(vn*(1:extension_distance));
uuc = uu+repmat(endpt.c,1,extension_distance);
vvr = vv+repmat(endpt.r,1,extension_distance);

labeled_endpts = threshed.lm;
labeled_endpts(~endpts) = 0;
neighborsFilter = ones(3);
neighborsFilter(2,2) = 0;
extended = double(threshed.lm > 0);
keepExtending = ~isnan(u);
extend_idx = sub2ind([1024 1024],vvr,uuc);
extended_diag = bwmorph(threshed.lm > 0,'diag',Inf);
% keepExtending(keepExtending) = any(extended_diag(extend_idx(keepExtending,:)),2);
extend_neighbors = zeros(size(keepExtending));



% neighborsFilter(1,1) = 0;
% neighborsFilter(3,1) = 0;
% neighborsFilter(3,3) = 0;
% neighborsFilter(1,3) = 0;

for i=1:extension_distance
    neighbors = imfilter(extended,neighborsFilter);
%     neighbors = imfilter(threshed.lm > 0,eight);
    extend_idx = sub2ind([1024 1024],vvr(keepExtending,i),uuc(keepExtending,i));
    extended(extend_idx) = 1;
    extend_neighbors = zeros(size(keepExtending));
    extend_neighbors(keepExtending) = neighbors(extend_idx);
    keepExtending = keepExtending & extend_neighbors < 3;
end
extended = bwmorph(extended,'skel',Inf);
extended = bwmorph(extended,'spur',Inf);
figure;
imshowpair(extended,threshed.lm > 0);
figure;
regions = bwlabel(~extended,4);
regions(regions == 1) = 0;
showPropMatrix(regions);

%% Figures
Ia = imadjust(I,[0.2 0.8]);
Ihat = imsubtract(imadd(Ia,imtophat(Ia,strel('disk',5))),imbothat(Ia,strel('disk',5)));
Ihat_skel = shiftdim (Ihat_skel,2);
Ihat_skel(1,M) = 1;
Ihat_skel(2,M) = 0;
Ihat_skel(3,M) = 0;

