function [trackModels res] = getTrackModels(allTrackParams, nPoints)

nTracks = numel(nPoints);

trackModels = zeros(nTracks, 4);
res = cell(nTracks, 1);

% pFirst and pLast are indexing allTrackParams
pLast = cumsum(nPoints);
pFirst = pLast-nPoints+1;

for iTrack = 1:nTracks
    params = allTrackParams(pFirst(iTrack):pLast(iTrack), [1 2 4 6]);
    params = num2cell(params,1);
    [x, y, sx, t] = params{:};
    ct = cos(t);
    st = sin(t);
    x = [x + sx .* ct; x - sx .* ct];
    y = [y + sx .* st; y - sx .* st];
    
    [trackModels(iTrack, :), res{iTrack}] = getSegmentModel(x,y);
end
