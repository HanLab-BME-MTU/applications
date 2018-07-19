function savePairTracks(movieData, allFeatures, tFirst, E, pFirst1, pFirst2, filename)

nFrames = movieData.nImages(1);
tOverlapFirst = max(tFirst(E(:,1)), tFirst(E(:,2)));

pairTracks = cell(nFrames,1);

for iFrame = 1:nFrames
    isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
    
    % Get the coordinates of the extremities of the valid pair
    offset = iFrame - tOverlapFirst(isPairInFrame);
    ind1 = pFirst1(isPairInFrame) + offset;
    ind2 = pFirst2(isPairInFrame) + offset;
    
    x1 = allFeatures(ind1,1);
    x2 = allFeatures(ind2,1);
    y1 = allFeatures(ind1,2);
    y2 = allFeatures(ind2,2);
    
    pairTracks{iFrame} = [x1 y1 x2 y2];
end

save(filename, 'pairTracks');
