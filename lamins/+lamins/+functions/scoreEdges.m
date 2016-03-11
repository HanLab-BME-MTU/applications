function [ score ] = scoreEdges( S,I )
%scoreEdges Scores the edges using distance weighted intensity metrics

import lamins.functions.*;

[~, EF] = S.faceEdges;
numFacesPerEdge = cellfun('length',EF);
edgeFilter = numFacesPerEdge == 2;
twoFaceEdges = find(edgeFilter);
EE = [EF{edgeFilter}]';
sets = independentSetDecomposition(EE);
score = zeros(1,S.edges.NumObjects);
out = cell(1,length(sets));
better = cell(1,length(sets));
for s=1:length(sets)
    better{s} = twoFaceEdges(sets{s});
end
better{end+1} = find(numFacesPerEdge == 1);
for s=1:length(better)
%     s
%     score(twoFaceEdges(sets{s})) = deleteEdgesAndScore(S,I, twoFaceEdges(sets{s}));

    out{s} = deleteEdgesAndScore(S,I, better{s});
end
for s=1:length(better)
    score(better{s}) = out{s};
end
% NaN?
% score(~edgeFilter) = NaN;

end
