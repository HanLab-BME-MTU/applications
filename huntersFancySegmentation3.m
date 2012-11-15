function mask = huntersFancySegmentation3(imageIn,varargin)
%THE FINAL CHAPTER!!! THE BEST OF BOTH WORLDS!! PREPARE TO BE AMAZED!!!

ip = inputParser;
ip.addParamValue('SigmasXY',[1 2 4],@(x)(all(x>=1)));
ip.addParamValue('SigmasZ',[1 2 4],@(x)(all(x>=1)));
ip.addParamValue('WeightZ',2,@(x)(numel(x) == 1 && x > 0));
ip.parse(varargin{:});
p = ip.Results;



nSTDintensity = 3;
%nSTDsurface = 2.5;%Yeah, like 2 or 3 are any less arbitrary because they're integers!! Fuck you, it works!
nSTDsurface = 2;

showImarisPlots = false;

maxResp = multiscaleSurfaceFilter3D(imageIn,p);

[backMean,backSTD] = robustMean(double(imageIn(:)),[],2);

%Use these background statistics to create a background mask
backThresh = backMean + nSTDintensity*backSTD;
backMask = imageIn > backThresh;

%Now get the high foreground threshold and foreground mask. This is high enough to remove some
%dim branches, but also excludes scattered light outside the cell.
foreThresh = thresholdOtsu(imageIn(backMask(:)));

foreMask = imageIn>foreThresh;

surfBackMean = mean(maxResp(:));
surfBackSTD = std(maxResp(:));

surfThresh = surfBackMean + (nSTDsurface*surfBackSTD);



surfMask = maxResp > surfThresh;


mask = surfMask | foreMask;


if showImarisPlots
   
    imMax = max(imageIn(:));
    imarisShowArray(cat(5,imageIn,uint16(foreMask)*imMax,uint16(surfMask)*imMax,uint16(mask)*imMax));
    
end

