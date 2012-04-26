function bp = analyze3DImageMaskedIntensities(images,mask,skelGraph,maskProp,curvSampRad,nucMask)
%ANALYZE3DIMAGEMASKEDINTENSITIES calculates several statistics regarding the spatial variation of the image intensities within the input mask 
% 
% bp = analyze3DImageMaskedIntensities(images,mask,skelGraph,maskProp)
% 
% This function analyzes how the intensities in the input image(s) vary
% within the input mask, including how they vary within each
% branch/skeletal element defined by the input skeletons.
% 
% Input:
% 
%   images - MxNxPxI matrix of concatenated 3D images, where the last
%   dimension corresponds to different images.
% 
%   mask - MxNxP logical matrix containing the mask.
% 
%   skelGraph - structure describing skeleton graph, as produced by
%   skel2graph.m
% 
%   maskProp - structure describing the mask's geometric properties, as
%   produced by analyze3DMaskGeometry.m
% 
%   maskPP - MxNxP logical matrix containing true values corresponding to
%   pixels in the mask which were added by post-processing.
%
%   curvSampRad - The radius, in pixels, to use when sampling the local
%   intensity and surface curvature.
%   Optional. Default is sqrt(2) which gives a 18-connected neighborhood.
%
% Output:
% 
%   bp - Structure containing the various calculated intensity statistics.
%        Most of the fields are self-explanatory, but I should probably add
%        a better description here someday...
%
%
% Hunter Elliott
% 1/2012 ... and I'm still in grad school...
%

%% -------------- Parameters ----------------- %%

distBinSz = 1;%Distance bin size for intensity vs. distance histograms.
distOut = 10;%Max distance outside the mask to calculate intensity stats
nIntBins = 100;%Number of bins for intensity histograms

%% ----------------- Input ------------------ %%

if nargin < 4
    error('Not enough inputs! Input some more stuff!')
end
    
if nargin < 5 || isempty(curvSampRad)
    curvSampRad = sqrt(2);%This gives the 18-connected neighborhood
end

if nargin < 6 || isempty(nucMask)
    nucMask = false(size(mask));
end

[M,N,P,nChan] = size(images);

%TEMP - the rest of this could be optimized for speed some more...!!

%--Get distance transform of mask and background and combine them--%

distXIn = bwdist(~mask);
distXOut = bwdist(mask);
distX = distXIn;%Combine, inverting so that distances outside the mask are negative
distX(~mask) = -distXOut(~mask);

%% -------------- Whole-Mask Curvature vs. Intensity Analysis ----------- %%

%First we create the depht-normalized forms of the images. Due to the
%sampling method, the variation with distance from the membrane will allow
%the local geometry to affect the intensity sampling. These allow us to
%separate the effects of depth variation from variation across the cell
%surface.
depthNormImages = zeros(size(images));
for j = 1:nChan
    
    %This seems really dumb but I guess its the most practical way to be
    %able to use the mask to index this...
    tmp = images(:,:,:,j);
    maskedMin = min(tmp(mask(:)));
    maskedMax = max(tmp(mask(:)));
    
    %Normalize the intensity with respect to depth
    tmp = depthNormalizeImage(tmp,distXIn,nucMask);
    %And then restore the original masked intensity ranges, to ease
    %comparison later on and allow us to use the same histogram bins
    tmp(mask(:)) = scaleContrast(tmp(mask(:)),[],double([maskedMin maskedMax]));
    
    depthNormImages(:,:,:,j) = tmp;
end

nSurfFaces = size(maskProp.SmoothedSurface.faces,1);

%Find surface pixels using the distance transform
surfPix = distX == 1;

%Convert these to XYZ coordinates to agree with the surface mesh
surfPixI = find(surfPix);
[surfPixXYZ(:,2),surfPixXYZ(:,1),surfPixXYZ(:,3)] = ind2sub([M N P],surfPixI);
nSurfPix = size(surfPixXYZ,1);
%We also need to get XYZ coordinates for all the pixels so we can
%associate by proximity to the surface
bodyPixI = find(mask & ~nucMask);
[bodyPixXYZ(:,2),bodyPixXYZ(:,1),bodyPixXYZ(:,3)] = ind2sub([M N P],bodyPixI);

%And find the face centroid locations since we have curvature data at the
%per-face level.
faceCenters = zeros(nSurfFaces,3);
for j = 1:nSurfFaces
    faceCenters(j,:) = mean(maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(j,:),:),1);
end

%Because of the way we smooth the mesh, the distance between the surface
%voxel and the surface mesh will vary. This will also vary in a way that is
%dependent on the local geometry. So, to prevent this from biasing our data
%we first find the closest surface face.
[closestSurfFaceIdxs,closestSurfFaceDist] = KDTreeClosestPoint(faceCenters,surfPixXYZ);
%Now we find the mesh surface faces within the sampling radius of this
%closest point
surfFaceIdxs = KDTreeBallQuery(faceCenters,faceCenters(closestSurfFaceIdxs,:),curvSampRad);
%For the image sampling, we use the actual surface voxel and find voxels
%within the specified distance. The local geometry will still affect the
%size of the sample (nPixels), but at least not the depth within the cell
%we are sampling, so this should still be unbiased
samplePixIdxs = KDTreeBallQuery(bodyPixXYZ,surfPixXYZ,curvSampRad);

bp.nPixPerCurvSamp = nan(nSurfPix,1);
bp.nFacesPerCurvSamp = nan(nSurfPix,1);
bp.dToSurfSamp = closestSurfFaceDist;
bp.surfPixInd = surfPixI;
bp.maxDepthPerCurvSamp = nan(nSurfPix,1);
bp.meanDepthPerCurvSamp = nan(nSurfPix,1);
bp.gaussCurvSampMean = nan(nSurfPix,1);
bp.gaussCurvSampSTD = nan(nSurfPix,1);
bp.meanCurvSampMean = nan(nSurfPix,1);
bp.meanCurvSampSTD = nan(nSurfPix,1);
bp.PC1CurvSampMean = nan(nSurfPix,1);
bp.PC1CurvSampSTD = nan(nSurfPix,1);
bp.PC2CurvSampMean = nan(nSurfPix,1);
bp.PC2CurvSampSTD = nan(nSurfPix,1);
bp.MaxAbsPCCurvSampMean = nan(nSurfPix,1);
bp.MaxAbsPCCurvSampSTD = nan(nSurfPix,1);
bp.intForCurvSampMean = nan(nSurfPix,nChan);
bp.intForCurvSampSTD = nan(nSurfPix,nChan);
bp.intForCurvSampMin = nan(nSurfPix,nChan);
bp.intForCurvSampMax = nan(nSurfPix,nChan);
bp.intForCurvSampMeanDepthNorm = nan(nSurfPix,nChan);
bp.intForCurvSampSTDDepthNorm = nan(nSurfPix,nChan);
bp.intForCurvSampMinDepthNorm = nan(nSurfPix,nChan);
bp.intForCurvSampMaxDepthNorm = nan(nSurfPix,nChan);


nMaskPix = numel(mask);%Needed later for indexing

%Calculate this ahead of time to avoid redundant calculations below
maxPCCurv = max(abs(real(maskProp.CurvaturePC1)),abs(real(maskProp.CurvaturePC2)));

%No go through the samples and calculate the averages etc.
for j = 1:nSurfPix                  
              
    %Get the indices of the current image sample voxels to save typing
    currInd = bodyPixI(samplePixIdxs{j});
    
    bp.nPixPerCurvSamp(j) = numel(currInd);
    bp.nFacesPerCurvSamp(j) = numel(surfFaceIdxs{j});
    
    bp.gaussCurvSampMean(j) = mean(maskProp.GaussianCurvature(surfFaceIdxs{j}));
    bp.gaussCurvSampSTD(j) = mean(maskProp.GaussianCurvature(surfFaceIdxs{j}));
    bp.meanCurvSampMean(j) = mean(maskProp.MeanCurvature(surfFaceIdxs{j}));
    bp.meanCurvSampSTD(j) = std(maskProp.MeanCurvature(surfFaceIdxs{j}));
    bp.PC1CurvSampMean(j) = mean(maskProp.CurvaturePC1(surfFaceIdxs{j}));
    bp.PC1CurvSampSTD(j) = std(maskProp.CurvaturePC1(surfFaceIdxs{j}));
    bp.PC2CurvSampMean(j) = mean(maskProp.CurvaturePC2(surfFaceIdxs{j}));
    bp.PC2CurvSampSTD(j) = std(maskProp.CurvaturePC2(surfFaceIdxs{j}));
    bp.MaxAbsPCCurvSampMean(j) = mean(maxPCCurv(surfFaceIdxs{j}));
    bp.MaxAbsPCCurvSampSTD(j) = std(maxPCCurv(surfFaceIdxs{j}));
    
    %Sample the fluorescence in each channel and the distance transform at these points                
    bp.maxDepthPerCurvSamp(j) = max(distX(currInd));        
    bp.meanDepthPerCurvSamp(j) = mean(distX(currInd));                
    for k = 1:nChan%Faster to do it this way, or to sample all channels at once and then re-shape the array?
        bp.intForCurvSampMean(j,k) = mean(double(images(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampSTD(j,k) = std(double(images(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampMin(j,k) = min(double(images(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampMax(j,k) = max(double(images(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampMeanDepthNorm(j,k) = mean(double(depthNormImages(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampSTDDepthNorm(j,k) = std(double(depthNormImages(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampMinDepthNorm(j,k) = min(double(depthNormImages(currInd + ((k-1)*nMaskPix))));
        bp.intForCurvSampMaxDepthNorm(j,k) = max(double(depthNormImages(currInd + ((k-1)*nMaskPix))));
    end      
end

%% ------------ Whole-Mask Intensity Vs Distance Analysis -------------- %%

%Set up the bins for the average intensity-vs-distance calculation
distBinEdges = -distOut:distBinSz:max(distX(:));
distBinCenters = distBinEdges(1:end-1) + distBinSz/2;
nDistBins = numel(distBinEdges)-1;

bp.wholeMaskMeanVsDist = nan(nChan,nDistBins);
bp.wholeMaskSTDVsDist = nan(nChan,nDistBins);
bp.wholeMaskNumPixVsDist = nan(nDistBins,1);
bp.wholeMaskDists = distBinCenters;
bp.wholeMaskIntHistVsDist = nan(nChan,nDistBins,nIntBins);%Histograms of intensity at each distance
bp.wholeMaskNormIntHistVsDist = nan(nChan,nDistBins,nIntBins);%Histograms of intensity at each distance, normalized at each distance to account for the varying # pixels there
bp.wholeMaskIntHistBins = nan(nChan,nIntBins);%Bins for these histograms

%TEMP - this can be sped up a bit - convert to index-arithmetic based
%sampling (and with branch sampling below!)
for j = 1:nDistBins                                            
    
    %Find the pixels for the current distance bin
    if distBinEdges(j) < 0
        currPix = distX >= distBinEdges(j) & distX < distBinEdges(j+1) & ~nucMask;
    else
        currPix = distX > distBinEdges(j) & distX <= distBinEdges(j+1) & ~nucMask;
    end
    bp.wholeMaskNumPixVsDist(j) = nnz(currPix);
    
    for k = 1:nChan
    
        tmp = double(images(:,:,:,k));%Lazy - makes indexing easier, decreases number of casts. Better way to the logical indexing below with diff dimension?  

        if j == 1
            %Set up intensity histogram bins for this channel
            bp.wholeMaskIntHistBins(k,:) = linspace(min(tmp(distX >= distBinEdges(1) & distX <= distBinEdges(end))),...
                                                 max(tmp(distX >= distBinEdges(1) & distX <= distBinEdges(end))),...
                                                 nIntBins);
        end
        
        %Get pixel stats for this distance bin        
        bp.wholeMaskMeanVsDist(k,j) = mean(tmp(currPix));
        bp.wholeMaskSTDVsDist(k,j) = std(tmp(currPix));        
        bp.wholeMaskHistVsDist(k,j,:) = histc(tmp(currPix(:)),bp.wholeMaskIntHistBins(k,:));
        bp.wholeMaskNormHistVsDist(k,j,:) = bp.wholeMaskHistVsDist(k,j,:) / bp.wholeMaskNumPixVsDist(j);

    end          
end


%% ---------------- Per-Branch Analysis ---------------- %%

[iTipVert,iTipEdge] = findTips(skelGraph.edges,size(skelGraph.vertices,1));
nTips = numel(iTipVert);

bp.branchTipPixelLenVsDepth = cell(nTips,1);%length vs. depth data for each pixel in each branch tip
bp.branchTipPixelInt = cell(nTips,1);%Intensity in each channel for each pixel in each branch tip
bp.branchTipPixelIntDepthNorm = cell(nTips,1);%Depth-normalized Intensity in each channel for each pixel in each branch tip

%Find the closest mask surface pixel to the branch tips, because some of these are
%non-integer valued, and if we just round them we end up outside the mask
%occasionally
surfIndForTips = KDTreeClosestPoint(surfPixXYZ,skelGraph.vertices(iTipVert,[2 1 3]));
surfIndForTips = surfPixI(surfIndForTips);

for j = 1:nTips
                
    %Get the geodesic distance in the mask from this branch tip
    currGDist = bwdistgeodesic(mask,surfIndForTips(j),'quasi-euclidean');
    %And find the geodesic distance to its base
    iBaseVert = skelGraph.edges(iTipEdge(j),:);
    iBaseVert = iBaseVert(iBaseVert ~= iTipVert(j));
    baseSub = round(skelGraph.vertices(iBaseVert,:));
    baseGeoDist = currGDist(baseSub(1),baseSub(2),baseSub(3));
    %Find pixels associated with this branch tip using this base distance.
    %This is only approximate and may include some non-branch areas of the mask.
    currBranchPix = currGDist <= baseGeoDist & ~nucMask;        
    %Extract the intensity, depth and length data for this branch in each
    %channel
    bp.branchTipPixelLenVsDepth{j} = nan(nnz(currBranchPix),2);
    bp.branchTipPixelInt{j} = nan(nnz(currBranchPix),nChan);
    
    bp.branchTipPixelLenVsDepth{j}(:,1) = currGDist(currBranchPix(:));%Length data for each pixel
    bp.branchTipPixelLenVsDepth{j}(:,2) = distX(currBranchPix(:));%Depth data for each pixel           
    
    for k = 1:nChan                
        tmpIm = images(:,:,:,k);%Seriously, there has to be a better way to do the logical indexing given the differing dimensionality here?? Right??
        bp.branchTipPixelInt{j}(:,k) = tmpIm(currBranchPix(:));
        tmpIm = depthNormImages(:,:,:,k);
        bp.branchTipPixelIntDepthNorm{j}(:,k) = tmpIm(currBranchPix(:));
    end                
    
%     %TMP - for testing and debug
%     imarisShowArray(cat(5,images(:,:,:,1),images(:,:,:,2),images(:,:,:,3),double(currGDist),double(currBranchPix)*100));    
    
end


    
