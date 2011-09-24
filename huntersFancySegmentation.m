function mask = huntersFancySegmentation(imageIn,maxResp)
% 
%Under construction.........
% 
%This one works okay for a decent number of images/movies, but requires
%this arbitrary enhancement constant, and leaves the cell-center empty
%which is difficult (impossible?) to completely fill in post-processing
%Also, it's not clear if


% sFiltSigma = 1;
% gFiltSigma = 1;
nSTDaboveBack = 3;
binSz = 10;
surfBinSz = 10;
showImarisPlots = false;
showPlots = false;
postProcess = false;

% %Background selection. The vast majority of the voxels are background, so
% %we approximate background values using the whole image.
% THis works alright, but tends to error on the high side
% backMean = mean(double(imageIn(:)));
% backSTD = std(double(imageIn(:)));

%Use histogram peak and lower width-at-half max
%This does well with mean but tends to under-estimate the STD
bins = double(min(imageIn(:)):binSz:max(imageIn(:)));
cts = histc(double(imageIn(:)),bins);
% [maxCt,iMax] = max(cts); %Approximate background mean as most common pixel value
% backMean = bins(iMax);
% [~,iHalfMax] = min(abs(maxCt/2 - cts(1:iMax)));%Assume that at least the lower half of the distribution is gaussian. We find the FWHM and use this to approximate
% backSTD = 2*(backMean-bins(iHalfMax)) / (2*sqrt(2*log(2)));%Get approximate sigma from FWHM
 
[backMean,backSTD] = robustMean(double(imageIn(:)),[],2);

%Use these background statistics to creat a background mask
backThresh = backMean + nSTDaboveBack*backSTD;
backMask = imageIn > backThresh;

backMaskDilated = imdilate(backMask,binarySphere(4));
sampleBackCts = histc(double(imageIn(~backMaskDilated(:))),bins);
sampleBackSTD = std(double(imageIn(~backMaskDilated(:))));
sampleBackMean = mean(double(imageIn(~backMaskDilated(:))));

%Now get the high foreground threshold. This is high enough to remove some
%dim branches, but also excludes scattered light outside the cell.
foreThresh = thresholdOtsu(imageIn(backMask(:)));

foreMask = imageIn>foreThresh;

%Use k-means to cluser pixels based on 
% [clustInd,centRoids] = kmeans(double(imageIn(backMask(:))),2);
% foreMask = zeros(size(imageIn));
% foreMask(backMask(:)) = clustInd;
% foreMask = foreMask == 2;

%[maxResp,d2X,d2Y,d2Z,maxRespScale] = multiscaleSurfaceFilter3D(imageIn);

% %ENHANCE!!! ENHANCE!! - This works, but In the end there needs to be a
% %better way to combine these two images........(i.e. less arbitrary!!)
enhancementFactor = 3*pi;% ;)
enhancedImage = double(imageIn) + maxResp * enhancementFactor;
surfThresh = 11; %All the way to 11


%Print the damn thing!
mask = enhancedImage > foreThresh;

%surfThresh = thresholdOtsu(maxResp(backMask(:))); %Doesn't work - selects too high
%surfThresh = thresholdOtsu(maxResp(foreMask(:))); %Doesn't work - selects too high

% [surfMean,surfSTD] = robustMean(maxResp(backMask(:)),[],1.8);
% surfThresh = surfMean + 3*surfSTD; %Works sometimes, other times too low.



%mask = imageIn > foreThresh | maxResp > surfThresh & backMask;

if showPlots
        
    disp(['Back Thresh: ' num2str(backThresh) ', Fore Thresh: ' num2str(foreThresh) ', Surf Thresh: ' num2str(surfThresh)])
    
    figure;    
    loglog(bins,cts);        
    hold on
    plot(ones(1,2)*backMean,ylim,'--r');
    plot(ones(1,2)*(backMean+backSTD),ylim,'--r');
    plot(ones(1,2)*backThresh,ylim,'--m');            
    
    plot(bins,sampleBackCts*(max(cts)/max(sampleBackCts)),'g');
    plot(ones(1,2)*(sampleBackMean),ylim,'--g')
    plot(ones(1,2)*(sampleBackMean+sampleBackSTD),ylim,'--g')
    
    plot(ones(1,2)*foreThresh,ylim,'--k')
        
    xLim = xlim;
    yLim = ylim;
    yNorm = max(cts)*exp(-(bins-backMean) .^2 ./ (2*backSTD^2));
    plot(bins,yNorm,'r')
    xlim(xLim);
    ylim(yLim);
    
    surfBins = min(maxResp(:)):surfBinSz:max(maxResp(:));
    surfCts = histc(maxResp(:),surfBins);
    surfBackCts = histc(maxResp(~backMask(:)),surfBins);
    surfForeCts = histc(maxResp(backMask(:)),surfBins);
    surfHighCts = histc(maxResp(foreMask(:)),surfBins);
    
    figure
    loglog(surfBins,surfCts ./ max(surfCts))
    hold on;
    plot(surfBins,surfBackCts ./ max(surfBackCts),'g');
    plot(surfBins,surfForeCts ./ max(surfForeCts),'r');
    plot(surfBins,surfHighCts ./ max(surfHighCts),'m');
    plot(ones(1,2)*surfThresh,ylim,'--b')
    if ~showImarisPlots
        figure;
        spy3d(backMask);
    end
end

maskUnProc = mask;

mask = postProcess3DMask(mask);

if showImarisPlots
    imarisShowArray(cat(5,imageIn,double(backMask)*max(imageIn(:)),double(foreMask)*max(imageIn(:)),double(maskUnProc)*max(imageIn(:)),double(mask)*max(imageIn(:))))    
end

%Gradient filtering
% [dX,dY,dZ] = gradientFilterGauss3D(double(imageIn),gFiltSigma);
% gMag = sqrt(dX .^2 + dY .^2 + dZ .^2);
% gradNMS = nonMaximumSuppression3D(dX,dY,dZ);
% 
% gradNMS(~backMask) = 0; %Remove background gradient values
% gnmsCC = bwconncomp(gradNMS>0); %Label the remaining components
% 
% gradPixPerCC = cellfun(@numel,gnmsCC.PixelIdxList);%Number of pixels per object
% [~,gradLargeToSmallCC] = sort(gradPixPerCC,'descend');
% gradSelectedCC = false(size(imageIn));
% gradSelectedCC(gnmsCC.PixelIdxList{gradLargeToSmallCC(1)}) = true;
% 
% jkl=1;

% if showPlots
%     bins = min(gMag(gradSelectedCC(:))):5:max(gMag(gradSelectedCC(:)));
%     cts = histc(gMag(gradSelectedCC(:)),bins);
%     figure
%     plot(bins,cts ./ max(cts));
%     hold on
%     bins = min(gMag(intMask(:))):5:max(gMag(intMask(:)));
%     cts = histc(gMag(intMask(:)),bins);
%     plot(bins,cts ./ max(cts),'r')
% end
% 
% gradThresh = min(gMag(gradSelectedCC(:)));
% 
% % %Surface filtering - enhances filaments and membranes
% % [d2X,d2Y,d2Z] = surfaceFilterGauss3D(imageIn,sFiltSigma);
% % %We only want areas with negative curvature, so we remove positive values
% % %prior to nms
% % d2X(d2X>0) = 0;
% % d2Y(d2Y>0) = 0;
% % d2Z(d2Z>0) = 0;
% % sMag = sqrt(d2X .^2 + d2Y .^2 + d2Z .^2);%Magnitude of surface filter response
% % %Apply non-maximum suppression and post-process.
% % surfNMS = nonMaximumSuppression3D(d2X,d2Y,d2Z);
% % %Remove background values
% % surfNMS(~intMask) = 0;
% % surfNMS(gradNMS>0) = 0;
% % 
% % %Post-process the NMS to select areas associated with the cell membrane.
% % snmsCC = bwconncomp(surfNMS>0);
% % surfPixPerCC = cellfun(@numel,snmsCC.PixelIdxList);%Number of pixels per object
% % [~,surfLargeToSmallCC] = sort(surfPixPerCC,'descend');
% % surfSelectedCC = false(size(imageIn));
% % surfSelectedCC(snmsCC.PixelIdxList{surfLargeToSmallCC(1)}) = true;
% 
% 
% %intMask(gradSelectedCC) = false;
% 
% mask = gMag > gradThresh;
% 
% if showImarisPlots
%     
%     imarisShowArray(cat(5,imageIn,gMag,double(mask)*1e3,double(intMask)*1e3));
%     
% end
% if showPlots
%     figure
%     spy3d(gradSelectedCC);    
%     figure
%     p = isosurface(gMag,gradThresh);
%     patch(p,'FaceAlpha',.4,'EdgeColor','none');
%     axis image,light
% end
% 



