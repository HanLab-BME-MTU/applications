function score = deleteEdgesAndScore(S,I,e)
% deleteEdgesAndScore Delete edges and score them
    import connectedComponents.*;
    import lamins.functions.*;
    S.faces = [];
    S2 = S.copy;
    S2.deleteEdges(e);
%     S2.bw = [];
%     S2.faces = [];
    faceStack = cat(3,S.faces.lm,S2.faces.lm);
    uFaceStack = unique(reshape(faceStack,[],2),'rows');
    uFaceStack = uFaceStack(uFaceStack(:,1) & uFaceStack(:,2),:);
%     F1 = S.getFaceProperties(I);
%     F2 = S2.getFaceProperties(I);
    F1 = S.getDistanceWeightedIntensity(I);
    F2 = S2.getDistanceWeightedIntensity(I);
%             distanceWeightedIntensity = cellfun(@(x) sum(D(x).*I(x))./sum(D(x)), faces.PixelIdxList , 'UniformOutput', false);

    
    cost = [F2.DistanceWeightedIntensity]' - accumarray(uFaceStack(:,2),[F1(uFaceStack(:,1)).DistanceWeightedIntensity],[],@min);
    midpoints = getEdgeMidpoints(ccFilter(S.edges,e));
    newFaces = S2.faces.lm(midpoints);
    score = zeros(size(e));
    score(newFaces > 0) = cost(newFaces(newFaces > 0));
    % maybe NaN
    score(newFaces == 0) = NaN;
end