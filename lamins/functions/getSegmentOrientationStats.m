function [ out ] = getSegmentOrientationStats( theta, skel )
%getSegmentOrientationStats get the orientation of the segments

% get bw image
if(~islogical(skel))
    skel = skel ~= 0;
end
cc = bwconncomp(skel,8);
rp = regionprops(cc,theta,'MeanIntensity','Extrema','PixelValues','Area');

out.cc = cc;
out.rp = rp;

end

