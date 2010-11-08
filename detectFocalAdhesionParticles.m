function [x y amp] = detectFocalAdhesionParticles(ima, mask, sigmaPSF, kSigma)

ima = double(ima);
[nrows ncols] = size(ima);

% Filter image with laplacian
bandPassIso = filterLoG(ima,sigmaPSF);
bandPassIso(~mask) = 0;

% Filter image with ridge detector
[~, ~, nms] = steerableFiltering(ima,2,sigmaPSF);
nms(~mask) = 0;

% Filter image with a Gaussian filter with support to compute bkg
% values.
lowPass = Gauss2D(ima, sigmaPSF);

locMaxIso = locmax2d(bandPassIso, [5 5]);
locMaxAni = locmax2d(nms, [5 5]);

% Sample the image boundary

% Be careful on the border condition !!! Eventhough sampled point are on
% the image boundary, they needs to be locmin !!!!!!

ptsBoundary = [...
    (1:nrows)', repmat(1,nrows,1); ...
    repmat(nrows,ncols-2,1), (2:ncols-1)'; ...
    (nrows:-1:1)', repmat(ncols,nrows,1); ...
    repmat(1,ncols-2,1), (ncols-1:-1:2)'];

ind = sub2ind(size(ima),ptsBoundary(1:10:end,1),ptsBoundary(1:10:end,2));
[y x] = ind2sub(size(ima),ind(mask(ind) == true));
ptsBoundary = [x y];

% Sample the cell outline
bwDist = bwdist(1 - mask);
outline = contourc(double(bwDist), [0, 0]);
nChunks = 0;
pos = 1;
while pos < size(outline,2)
    n = outline(2,pos);
    pos = pos + n + 1;
    nChunks = nChunks + 1;
end
ptsCellEdge = cell(nChunks,1);
pos = 1;
for iChunk = 1:nChunks
    n = outline(2,pos);
    pts = outline(:,pos+1:pos+n);
    l = sum(sqrt(sum(diff(pts,1,2).^2,1)));
    s = spline(linspace(0,l,n), pts);
    ptsCellEdge{iChunk} = fnval(s,0:10:l)';
    pos = pos + n + 1;
end
ptsCellEdge = vertcat(ptsCellEdge{:});
    
% Get local minima from the low pass filtering image and the edge outline
% The definition of local minima is not good since a local maxima could
% happen right in the middle of a bundle of adhesion and therefore have a
% very high intensity value. We therefore need to select points that look
% uniform in a small patch area.

locMin = locmin2d(lowPass, [3 3]);
indMin = find(locMin & mask);
[yMin xMin] = ind2sub(size(locMin), indMin);
ptsMin = [xMin yMin];
ptsMin = vertcat(ptsBoundary, ptsCellEdge, ptsMin);

% Get local maxima
indMax = find(locMaxIso ~= 0 | locMaxAni ~= 0);
[qy qx] = ind2sub(size(ima), indMax);

% Compute expected background value at locMax position based on loc min
% values.
x = ptsMin(:,1);
y = ptsMin(:,2);
z = interp2(lowPass,ptsMin(:,1),ptsMin(:,2));

F = TriScatteredInterp(x,y,z);

[X Y] = meshgrid(floor(min(x)):ceil(max(x)), floor(min(y)):ceil(max(y)));
Z = F(X,Y);
mesh(X,Y,Z);
qz = F(qx,qy);
hold on;
plot3(qx,qy,qz,'o');

% Compute triangulation
dt = DelaunayTri(ptsMin(:,1),ptsMin(:,2));
iTri = pointLocation(dt,[qx, qy]);
assert(all(~isnan(iTri)));
tri = dt(iTri,:);

% For each point, compute the intensity on the triangle by bilinear
% interpolation

figure, imshow(ima,[]);
hold on;
plot(qx, qy, 'r.');
triplot(dt);

Delaunay

% 
% % Compute the adaptative mean
avgBkg = inf(size(indMax));
% avgBkg(nnzIdx) = mean(reshape(lowPass(indMin(tri(:))), size(tri)),2);
% 
% % Compute the adaptive standard deviation (over 27 points)
% stdBkg = inf(size(indMax));
% 
% ind3x3 = bsxfun(@plus,[-nrows 0 nrows],(-1:1)');
% ind3x3b1 = bsxfun(@plus,indMin(tri(:,1)),ind3x3(:)');
% ind3x3b2 = bsxfun(@plus,indMin(tri(:,2)),ind3x3(:)');
% ind3x3b3 = bsxfun(@plus,indMin(tri(:,3)),ind3x3(:)');
% ind3x3b = [ind3x3b1, ind3x3b2, ind3x3b3];
% 
% stdBkg(nnzIdx) = std(ima(ind3x3b),[],2);
% 
% finalIndMax = indMax(ima(indMax) > avgBkg + kSigma * stdBkg);
% 
% [y x] = ind2sub(size(ima), finalIndMax);
% 
% amp = ima(finalIndMax);
% 
% % Detection of particles on ridges