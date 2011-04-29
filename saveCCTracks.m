function saveCCTracks(movieData, CC, allFeatures, tFirst, tLast, pFirst, filename)

nFrames = movieData.nImages(1);
nCC = numel(CC);

% Compute the first and last frame of each CC
tFirstCC = cellfun(@(trackIdx) min(tFirst(trackIdx)), CC); % first frame of CC
tLastCC = cellfun(@(trackIdx) max(tLast(trackIdx)), CC);   % last frame of CC

colors = hsv(nCC);

segments = cell(nFrames,1);

for iFrame = 1:nFrames    
    isInFrame = find(iFrame >= tFirstCC & iFrame <= tLastCC);
    
    pFirstCC = cellfun(@(trackIdx) ...
        pFirst(trackIdx) + iFrame - tFirst(trackIdx), ...
        CC(isInFrame), 'UniformOutput', false);
    
    isTrackInFrame = cellfun(@(trackIdx) iFrame >= tFirst(trackIdx) & ...
        iFrame <= tLast(trackIdx), CC(isInFrame), 'UniformOutput', false);
    
    % Gather every feature of each track in the CC between at iFrame
    allFeaturesCC = cellfun(@(aa,bb) arrayfun(@(a,b) ...
        allFeatures(a:a+b-1, [1 2 4 6]), aa, bb, 'UniformOutput', false), ...
        pFirstCC, isTrackInFrame, 'UniformOutput', false);
    
    allFeaturesCC = cellfun(@(c) vertcat(c{:}), allFeaturesCC, ...
        'UniformOutput', false);
    
    modelsCC = getSegmentModels(allFeaturesCC);
    
    x = modelsCC(:,1);
    y = modelsCC(:,2);
    l = modelsCC(:,3);
    t = modelsCC(:,4);
    ct = cos(t);
    st = sin(t);
    
    x1 = x + ct .* l/2;
    x2 = x - ct .* l/2;
    y1 = y + st .* l/2;
    y2 = y - st .* l/2;
    
    segments{iFrame} = [x1 y1 x2 y2 colors(isInFrame,:)];
end

save(filename, 'segments');
