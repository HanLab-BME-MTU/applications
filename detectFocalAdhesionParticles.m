function [x y amp] = detectFocalAdhesionParticles(ima, mask, sigmaPSF, kSigma)

ima = double(ima);
nrows = size(ima,1);

% Filter image with laplacian
bandPassIso = filterLoG(ima,sigmaPSF);
bandPassIso(~mask) = 0;

% Filter image with ridge detector
[bandPassAni, theta, nms] = steerableFiltering(ima,2,sigmaPSF);
bandPassAni(~mask) = 0;
nms(~mask) = 0;

% Filter image with a Gaussian filter with large support to compute bkg
% values.
lowPass = Gauss2D(ima,4 * sigmaPSF);
lowPass(~mask) = 0;

% test
locMaxIso = locmax2d(bandPassIso, [5 5]);
locMaxAni = locmax2d(nms, [5 5]);
ind = find(locMaxIso ~= 0 | locMaxAni ~= 0);
[y x] = ind2sub(size(ima), ind);
amp = ima(ind);

% imshow(ima,[]); hold on;
% [y x] = ind2sub(size(ima), indIso);
% plot(x,y,'g.');
% [y x] = ind2sub(size(ima), indAni);
% plot(x,y,'b.');
% [y x] = ind2sub(size(ima), indBoth);
% plot(x,y,'r.');



% % Detection of isolated particles
% hside1 = 2 * ceil(3 * sigmaPSF) + 1;
% hside2 = 2 * ceil(sigmaPSF) + 1;
% locMax = locmax2d(bandPassIso, [hside1 hside1]);
% % Remove local maxima lying on ridges since these locations will be
% % evaluated later anyways.
% locMax(nms ~= 0) = 0;
% % Get local minima from the low pass filtering image
% locMin = locmin2d(lowPass, [hside2 hside2]);
% 
% indMin = find(locMin);
% indMax = find(locMax);
% [yMin xMin] = ind2sub(size(locMin), indMin);
% [yMax xMax] = ind2sub(size(locMax), indMax);
% dt = DelaunayTri(xMin,yMin);
% triangles = pointLocation(dt,[xMax, yMax]);
% nnzIdx = ~isnan(triangles);
% tri = dt(triangles(nnzIdx),:);
% 
% % Compute the adaptative mean
% avgBkg = inf(size(indMax));
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