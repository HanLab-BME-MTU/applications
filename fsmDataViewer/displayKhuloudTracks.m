function displayKhuloudTracks(hAxes, tag, tracks, layerColor)

nTracks = numel(tracks);

for iTrack = 1:nTracks
    pts = tracks{iTrack};
    line(pts(:,1), pts(:,2), 'Marker', 'none', 'Color', layerColor, 'Parent', hAxes, 'Tag', tag);
    line(pts(end,1),pts(end,2),'Marker','.','Color', layerColor, 'Parent', hAxes, 'Tag', tag);
end