function [skel] = laminSkeleton(I,mask)

    if(nargin < 2)
        mask = ones(size(I));
    end

I = im2double(imadjust(I,stretchlim(I,0)));

%% Do steerable detection
[res, theta, nms] = steerableDetector(I,4,5);
skel = nms ~= 0 & mask;
skel = bwmorph(skel,'skel',Inf);

%% Do thresholding by area
theta_stats = getSegmentOrientationStats(theta,skel);
areaThresh = thresholdRosin([theta_stats.rp.Area]);
areaThresh = 5;
threshed.cc = filtercc(theta_stats.cc,[theta_stats.rp.Area] > areaThresh);
threshed.lm = labelmatrix(threshed.cc);
skel = threshed.lm > 0;

%% Identify end points and branch points
endpts = bwmorph(skel,'endpoints');


vector = getEndPointVector(skel,'local');

maxPixelExtension = 20;
connectivity = 8;
skel = extendVectorUntilConnected(skel,endpts,vector,maxPixelExtension,connectivity);

if(nargout == 0)
    %imshow(skel,[]);
end

end
