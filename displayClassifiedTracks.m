function displayClassifiedTracks(hAxes, tag, tracks, layerColor)

nTracks = numel(tracks);

for iTrack = 1:nTracks
    pts = tracks(iTrack).trackCoords;
    color = tracks(iTrack).color;
    line(pts(:,1), pts(:,2), 'Marker', 'none', 'Color', color, 'Parent', hAxes, 'Tag', tag);
    line(pts(end,1),pts(end,2), 'MarkerSize', 10, 'Marker', '.', 'Color', 'r', 'Parent', hAxes, 'Tag', tag);
end