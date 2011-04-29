function displayClassifiedSegments(hAxes, tag, segments, layerColor)

nSegments = size(segments,1);

for iSeg = 1:nSegments
    
    pts = segments(iSeg,1:4);
    color = segments(iSeg,5:end);
    
    line(pts([1 3]), pts([2 4]), 'Marker', 'none', 'Color', color, 'Parent', hAxes, 'Tag', tag);
end