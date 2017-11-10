I = Lamins(2).cellReader{3,1,15};
I = im2double(imadjust(I,stretchlim(I,0)));
[res, theta, nms] = steerableDetector(double(i),4,5);
skel = nms ~= 0 & mask;
skel2 = bwmorph(imdilate(skel,strel('disk',2)),'skel',Inf);
antiskel = bwmorph(~imdilate(skel2,strel('disk',2)),'thin',Inf) & mask;