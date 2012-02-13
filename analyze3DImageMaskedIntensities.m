function bp = analyze3DImageMaskedIntensities(images,mask,skelGraph,maskProp,curvSampRad,nucMask,distX,distL,distXout,distOutL)
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
    nucMask = true(size(mask));
end

[M,N,P,nChan] = size(images);

%--Get distance transform of mask and background and combine them--%

% if numel(mask) < 2^24
%     [distX distL] = bwdist(~mask);%Internal mask distance transform (distance FROM background) and labelling by closest pixel
%     [distXOut distOutL] = bwdist(mask);%Background pixel distance transform - distance FROM mask    
% else    
%     %SHHHIIIITTTT!!! Due to a bug in bwdist which causes an incorrect label
%     %matrix to be returned when image has more than 2^24 elements (which in our
%     %data is invariably the case), we are forced to use the old version of
%     %bwdist which is SLOOWW. This is a known issue:
%     %http://www.mathworks.com/support/bugreports/737384 but currently there is
%     %not a fix and this is the workaround they suggest... !!
%     [distX distL] = bwdist_old(~mask);%Internal mask distance transform (distance FROM background) and labelling by closest pixel
%     [distXOut distOutL] = bwdist_old(mask);%Background pixel distance transform - distance FROM mask
%     distX(~mask) = -distXOut(~mask);%Combine, inverting so that distances outside the mask are negative
% end
% distX(~mask) = -distXOut(~mask);%Combine, inverting so that distances outside the mask are negative

%% -------------- Whole-Mask Curvature vs. Intensity Analysis ----------- %%


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
[surfFaceIdxs,surfFaceDist] = KDTreeBallQuery(faceCenters,faceCenters(closestSurfFaceIdxs,:),curvSampRad);
%For the image sampling, we use the actual surface voxel and find voxels
%within the specified distance. The local geometry will still affect the
%size of the sample (nPixels), but at least not the depth within the cell
%we are sampling, so this should still be unbiased
[samplePixIdxs,samplePixDist] = KDTreeBallQuery(bodyPixXYZ,surfPixXYZ,curvSampRad);

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
    end      
%     %TEMP - for debugging
%     if mod(j,10000) == 0
%         tmpMask = false(size(mask));
%         tmpMask(currInd) = true;        
%         imarisShowArray(cat(5,images(:,:,:,1),double(mask),double(tmpMask)))
%         figure
%         spy3d(mask,'.k'),hold on
%         plot3(faceCenters(surfFaceIdxs{j},1),faceCenters(surfFaceIdxs{j},2),faceCenters(surfFaceIdxs{j},3),'rx')
%         plot3(bodyPixXYZ(samplePixIdxs{j},1),bodyPixXYZ(samplePixIdxs{j},2),bodyPixXYZ(samplePixIdxs{j},3),'gx')
%         plot3(surfPixXYZ(j,1),surfPixXYZ(j,2),surfPixXYZ(j,3),'mo','MarkerSize',10)        
%     end
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

%TEMP - this can be sped up a bit...
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

for j = 1:nTips
            
    tipSub = round(skelGraph.vertices(iTipVert(j),:));%Get closes subscripts for non-integer tip location
    if ~mask(tipSub(1),tipSub(2),tipSub(3));
        %If not, use the background distance
        %labelling to find the closest mask surface pixel        
        [tipSub(1) tipSub(2) tipSub(3)] = ind2sub([M N P],distOutL(currPt));        
    end
    %Get the geodesic distance in the mask from this branch tip
    currGDist = bwdistgeodesic(mask,sub2ind([M N P],tipSub(1),tipSub(2),tipSub(3)),'quasi-euclidean');
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
    end                
    
end


    
