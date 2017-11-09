I = Lamins(2).cellReader{3,1,15};
I = im2double(imadjust(I,stretchlim(I,0)));
[res, theta, nms] = steerableDetector(double(I),4,5);
mask = Lamins(2).getNuclearMask(3,1,15);
skel = nms ~= 0 & mask;
skel = bwmorph(skel,'skel',Inf);
skel2 = bwmorph(imdilate(skel,strel('disk',2)),'skel',Inf);
antiskel = bwmorph(~imdilate(skel2,strel('disk',2)),'thin',Inf) & mask;