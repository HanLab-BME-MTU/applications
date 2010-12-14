function [x y amp] = detectFocalAdhesionParticles(ima, mask, sigmaPSF)

ima = double(ima);
%[nrows ncols] = size(ima);

% Filter image with laplacian
bandPassIso = filterLoG(ima,sigmaPSF);
bandPassIso(bandPassIso < 0) = 0;
bandPassIso(~mask) = 0;

% % Filter image with ridge detector
% bandPassAni = steerableFiltering(ima,2,sigmaPSF);
% bandPassAni(~mask) = 0;
% 
% % Compute the local maxima of the bandpass filtered images
locMaxIso = locmax2d(bandPassIso, [5 5]);

bw = blobSegmentThreshold(bandPassIso,0,0,mask);

locMaxIso(~bw) = 0;

indMax = find(locMaxIso);
[y x] = ind2sub(size(ima), indMax);

% locMaxAni = locmax2d(bandPassAni, [11 11]);
% 
% % Filter image with a Gaussian filter
% lowPass = Gauss2D(ima, sigmaPSF,1);
% %lowPass = padarray(lowPass,[1 1],'replicate');
% 
% % Reconstruct the cell boundary.
% mask = imdilate(mask,strel('square',3));
% bwDist = bwdist(~mask);
% outline = contourc(double(bwDist), [0, 0]);
% nChunks = 0;
% pos = 1;
% while pos < size(outline,2)
%     n = outline(2,pos);
%     pos = pos + n + 1;
%     nChunks = nChunks + 1;
% end
% ptsCellEdge = cell(nChunks,1);
% pos = 1;
% for iChunk = 1:nChunks
%     n = outline(2,pos);
%     ptsCellEdge{iChunk} = ceil(outline(:,pos+1:pos+n)');
%     pos = pos + n + 1;
% end
% ptsCellEdge = vertcat(ptsCellEdge{:});
% ptsCellEdge = ptsCellEdge( ...
%     ptsCellEdge(:,1) <= nrows & ...
%     ptsCellEdge(:,1) >= 1 & ...
%     ptsCellEdge(:,2) <= ncols & ...
%     ptsCellEdge(:,2) >= 1,:);
% 
% ind = sub2ind(size(ima),ptsCellEdge(:,2),ptsCellEdge(:,1));
% ptsCellEdge = ptsCellEdge(locmin1d(lowPass(ind),5),:);
% 
% % Get local minima from the low pass filtering image and the cell boundary
% lowPassSS = imresize(lowPass .* mask,.125);
% locMinSS = locmin2d(lowPassSS,[3 3]);
% indMinSS = find(locMinSS);
% [yMin xMin] = ind2sub(size(locMinSS), indMinSS);
% ptsMin = [xMin yMin] * 8;
% ptsMin = unique(vertcat(ptsCellEdge, ptsMin),'rows');
% indMin = sub2ind(size(lowPass),ptsMin(:,2),ptsMin(:,1));
% 
% % Get local maxima
% %indMax = find(locMaxIso ~= 0 | locMaxAni ~= 0);
% indMax = find(locMaxIso ~= 0);
% [yMax xMax] = ind2sub(size(ima), indMax);
% ptsMax = [xMax yMax];
% 
% % Compute the Delaunay triangulation
% dt = DelaunayTri(ptsMin(:,1),ptsMin(:,2));
% 
% % Compute which triangle enclosed local maxima
% iTri = pointLocation(dt,ptsMax);
% 
% isValid = ~isnan(iTri);
% indMax = indMax(isValid);
% ptsMax = ptsMax(isValid,:);
% iTri = iTri(isValid);
% 
% tri = dt(iTri,:);
% 
% % Compute the average value 
% avgBkg = mean(reshape(lowPass(indMin(tri(:))), size(tri)),2);
% 
% % Compute the adaptive standard deviation (over 27 points)
% ind3x3 = bsxfun(@plus,[-nrows 0 nrows],(-1:1)');
% ind3x3b1 = bsxfun(@plus,indMin(tri(:,1)),ind3x3(:)');
% ind3x3b2 = bsxfun(@plus,indMin(tri(:,2)),ind3x3(:)');
% ind3x3b3 = bsxfun(@plus,indMin(tri(:,3)),ind3x3(:)');
% 
% bkg1 = bsxfun(@minus,ima(ind3x3b1), mean(ima(ind3x3b1), 2));
% bkg2 = bsxfun(@minus,ima(ind3x3b2), mean(ima(ind3x3b2), 2));
% bkg3 = bsxfun(@minus,ima(ind3x3b3), mean(ima(ind3x3b3), 2));
% bkg = [bkg1, bkg2, bkg3];
% 
% stdBkg = std(bkg, [], 2);

% finalIndMax = indMax(ima(indMax) > avgBkg + kSigma * stdBkg);
% [y x] = ind2sub(size(ima), finalIndMax);
amp = ima(indMax);
