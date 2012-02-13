function nucMask = segment3DCellNucleus(image,mask)

%Completely ad-hoc method for segmenting out the nucleus based on location
%and intensity. This function sucks and should be made more general/less
%arbitrary!!!!!


showPlots = false;

safeZone = 4;%Buffer region at perimeter to exclude dark areas close to mask edge.
dRad = 4;
dilStrel = strel('ball',dRad,dRad,0);
dilIm = imdilate(image,dilStrel);

distX = bwdist(~mask);

centerMostPix = distX > prctile(distX(mask(:)),99);

intThresh = mean(image(centerMostPix(:))) + 2*std(image(centerMostPix(:)));


tmpMask = dilIm < intThresh & mask & distX > safeZone;

nmLabel = bwconncomp(tmpMask,6);
nPix = cellfun(@numel,nmLabel.PixelIdxList);
[~,iBySize] = sort(nPix,'descend');

nucMask = false(size(mask));
nucMask(nmLabel.PixelIdxList{iBySize(1)}) = true;

%Sort of half-assed reverse the effects of the dilation by dilating the
%mask
nucMask = imdilate(nucMask,strel('disk',dRad,0));

if showPlots
   
    imMax = max(image(:));
    imarisShowArray(cat(5,image,double(nucMask)*imMax/2,double(mask)*imMax/2))
    
    
end



