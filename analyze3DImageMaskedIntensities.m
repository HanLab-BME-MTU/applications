function bp = analyze3DImageMaskedIntensities(images,mask,skelGraph,maskProp)
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

%TEMP - right now this function is not really optimized for speed, and
%could be made quite a bit faster if necessary.... HLE

%Add some input checks? Who cares...

[M,N,P,nChan] = size(images);


%--Get distance transform of mask and background and combine them--%

%SHHHIIIITTTT!!! Due to a bug in bwdist which causes an incorrect label
%matrix to be returned when image has more than 2^24 elements, we are
%forced to use the old version of bwdist which is SLOOWW. This is a known
%issue: http://www.mathworks.com/support/bugreports/737384 but currently
%there is not a fix and this is the workaround they suggest... !!
% [distX distL] = bwdist(~mask);%Internal mask distance transform (distance FROM background) and labelling by closest pixel
% [distXOut distOutL] = bwdist(mask);%Background pixel distance transform - distance FROM mask
% distX(~mask) = -distXOut(~mask);%Combine, inverting so that distances outside the mask are negative

[distX distL] = bwdist_old(~mask);%Internal mask distance transform (distance FROM background) and labelling by closest pixel
[distXOut distOutL] = bwdist_old(mask);%Background pixel distance transform - distance FROM mask
distX(~mask) = -distXOut(~mask);%Combine, inverting so that distances outside the mask are negative



%% -------------- Whole-Mask Curvature vs. Intensity Analysis ----------- %%
nSurfPts = size(maskProp.SmoothedSurface.faces,1);
subSampRate = 50;%Takes WAAAAYYY too long to do every single surface point so we randomly sub-sample

rng('shuffle');%Init random number generator
iSamp = randperm(nSurfPts);
nSamp = floor(nSurfPts/subSampRate);
iSamp = iSamp(1:nSamp);

% iSamp = 1:subSampRate:nSurfPts;%TEMP - for debugging, just sample like this so we can reproduce errors etc.
% nSamp = numel(iSamp);

bp.intForCurvSamp = cell(nSamp,1);
bp.depthForCurvSamp = cell(nSamp,1);
bp.gaussCurvSamp = nan(nSamp,1);
bp.meanCurvSamp = nan(nSamp,1);

%Get the curvature data for these surface samples
bp.gaussCurvSamp = maskProp.GaussianCurvature(iSamp);
bp.meanCurvSamp = maskProp.MeanCurvature(iSamp);
    
nMaskPix = numel(mask);%Needed later for indexing

for j = 1:nSamp           
    
    %Find the closest mask surface pixel to this surface face (due to
    %smoothing, they aren't perfectly colocalized)        
    currPt = round(mean(maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(iSamp(j),:),:),1));%Use barycenter of face vertices
    %Convert to linear index to speed up the steps below
    currPt = sub2ind([M N P],currPt(2),currPt(1),currPt(3));
    %Due to the rounding above this point may not actually lie in the mask.
    if ~mask(currPt)
        %If not, use the background distance
        %labelling to find the closest mask surface pixel        
        currPt = distOutL(currPt);
    end        
    %TEMP - debug
    if ~mask(currPt)
        error('curr Pt not in mask1!!')
    end
    
    %Find pixels in the mask interior for which this is the closest surface
    %point using the distance labelling matrices
    
    %First, we find the closest background pixel to this mask pixel
    bakPt = distL(currPt);
    %Now we use this to find all mask pixels where this background pixel is
    %the closest background pixel.       
    closestPts = find(distL == distL(bakPt));%Logical indexing faster?? Even though it will be VERY sparse?????
    closestPts(closestPts == bakPt) = [];%Exclude the background pixel itself, which is found due to the indexing method.
    
    %closestPts = distL == distL(currPt(2),currPt(1),currPt(3));    
    %nPts = nnz(closestPts,1);        
    nPts = size(closestPts,1);
    
    if ~all(mask(closestPts))
        error('not all found points are in mask!')
    end            
    
    %Sample the fluorescence in each channel and the distance transform at these points                
    bp.depthForCurvSamp{j} = distX(closestPts);        
    bp.intForCurvSamp{j} = nan(nPts,nChan);
    for k = 1:nChan%Faster to do it this way, or to sample all channels at once and then re-shape the array?
        bp.intForCurvSamp{j}(:,k) = images(closestPts + ((k-1)*nMaskPix));
    end
    
%     %TEMP - for debugging
%     cp = false(size(mask));
%     cp(currPt) = true;
%     closePt = false(size(mask));
%     closePt(closestPts) = true;    
%     imMax = max(max(max(images(:,:,:,1))));
%     imarisShowArray(cat(5,double(mask)*imMax*2,double(cp)*imMax*2,double(closePt)*imMax*2,images(:,:,:,1)))
%     
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


    
