function mask = nuclearMaskFromNMS(nms)
% Uses the NMS from steerableDetector to create a nuclear mask
    bwNms = nms ~= 0;
    filled = imfill(bwNms,'holes') & ~bwNms;
    filled = imclose(filled,strel('disk',3));
    mask = filled;
end
