function bp = analyze3DImageMaskedIntensities(images,mask,skelGraph,maskProp,maskPP)
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

%TEMP - INPUT CHECKS!! DEFAULTS!!!

%TEMP - right now this function is not really optimized for speed, and
%could be made quite a bit faster if necessary.... HLE

%Add some input checks? Who cares...

[M,N,P,nChan] = size(images);

%--Get distance transform of mask and background and combine them--%

if numel(mask) < 2^24
    [distX distL] = bwdist(~mask);%Internal mask distance transform (distance FROM background) and labelling by closest pixel
    [distXOut distOutL] = bwdist(mask);%Background pixel distance transform - distance FROM mask    
else    
    %SHHHIIIITTTT!!! Due to a bug in bwdist which causes an incorrect label
    %matrix to be returned when image has more than 2^24 elements (which in our
    %data is invariably the case), we are forced to use the old version of
    %bwdist which is SLOOWW. This is a known issue:
    %http://www.mathworks.com/support/bugreports/737384 but currently there is
    %not a fix and this is the workaround they suggest... !!
    [distX distL] = bwdist_old(~mask);%Internal mask distance transform (distance FROM background) and labelling by closest pixel
    [distXOut distOutL] = bwdist_old(mask);%Background pixel distance transform - distance FROM mask
    distX(~mask) = -distXOut(~mask);%Combine, inverting so that distances outside the mask are negative
end
distX(~mask) = -distXOut(~mask);%Combine, inverting so that distances outside the mask are negative

%% -------------- Whole-Mask Curvature vs. Intensity Analysis ----------- %%

nSurfFaces = size(maskProp.SmoothedSurface.faces,1);

%Find surface pixels which are not the result of post-processing to be
%extra careful.
surfPix = distX == 1 & ~maskPP;

%Convert these to XYZ coordinates to agree with the surface mesh
surfPixI = find(surfPix);
[surfPixXYZ(:,2),surfPixXYZ(:,1),surfPixXYZ(:,3)] = ind2sub([M N P],surfPixI);
nSurfPix = size(surfPixXYZ,1);

%And find the face centroid locations since we have curvature data at the
%per-face level.
faceCenters = zeros(nSurfFaces,3);
for j = 1:nSurfFaces
    faceCenters(j,:) = mean(maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(j,:),:),1);
end

%Now find surface faces within the specified sampling radius of each of the
%surface pixels
[surfFaceIdxs,surfFaceDist] = KDTreeClosestPoint(faceCenters,surfPixXYZ);

%First, get the surface curv properties for this surface pixel
bp.gaussCurvSamp = maskProp.GaussianCurvature(surfFaceIdxs);
bp.meanCurvSamp = maskProp.MeanCurvature(surfFaceIdxs);
bp.PC1CurvSamp = maskProp.CurvaturePC1(surfFaceIdxs);
bp.PC2CurvSamp = maskProp.CurvaturePC2(surfFaceIdxs);  
bp.distToSurf = surfFaceDist;
bp.surfPixForCurvSamp = surfPixI;
bp.surfFaceForCurvSamp = surfFaceIdxs;

bp.intForCurvSamp = cell(nSurfPix,1);
bp.depthForCurvSamp = cell(nSurfPix,1);

nMaskPix = numel(mask);%Needed later for indexing
for j = 1:nSurfPix                  
            
    %tic          
    
    %Find pixels in the mask interior for which this is the closest surface
    %point using the distance labelling matrices
    
    %First, we find the closest background pixel to this mask pixel
    bakPt = distL(surfPixI(j));
    %Now we use this to find all mask pixels where this background pixel is
    %the closest background pixel.       
    closestPts = find(distL == distL(bakPt));%Logical indexing faster?? If we do it this way we can do the indexing trick below...
    %Now we want to exclude the background pixel itself, which is found due
    %to the indexing method, and any other surface pixels that were found -
    %these can be equidistant to this background point and so will be
    %returned by the label matrix, but they will have their own curvature
    %samples.
    closestPts(closestPts == bakPt |  distX(closestPts) == 1 & closestPts ~= surfPixI(j)) = [];
    
    nPts = size(closestPts,1);
        
    %Sample the fluorescence in each channel and the distance transform at these points                
    bp.depthForCurvSamp{j} = distX(closestPts);        
    bp.intForCurvSamp{j} = nan(nPts,nChan);
    for k = 1:nChan%Faster to do it this way, or to sample all channels at once and then re-shape the array?
        bp.intForCurvSamp{j}(:,k) = images(closestPts + ((k-1)*nMaskPix));
    end
        
%     %TEMP - for debugging
%     if nPts > Inf
%         cp = false(size(mask));
%         cp(surfPixXYZ(j,2),surfPixXYZ(j,1),surfPixXYZ(j,3)) = true;
%         closePt = false(size(mask));
%         closePt(closestPts) = true;    
%         imMax = max(max(max(images(:,:,:,1))));
%         imarisShowArray(cat(5,double(mask)*imMax*2,double(cp)*imMax*2,double(closePt)*imMax*2,images(:,:,:,3)))
%     end
%     
    if mod(j,round(nSurfPix/20))==0
        disp(j/nSurfPix)
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

%TEMP - this can be sped up a bit...
for j = 1:nDistBins                                            
    
    %Find the pixels for the current distance bin
    if distBinEdges(j) < 0
        currPix = distX >= distBinEdges(j) & distX < distBinEdges(j+1);
    else
        currPix = distX > distBinEdges(j) & distX <= distBinEdges(j+1);
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
bp.branchTipPixelInt = cell(nTips,1);%Insity in each channel for each pixel in each branch tip

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
    currBranchPix = currGDist <= baseGeoDist;        
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


    
