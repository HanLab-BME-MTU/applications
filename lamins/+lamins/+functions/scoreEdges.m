function [ score ] = scoreEdges( S,I )
%scoreEdges Scores the edges using distance weighted intensity metrics

import lamins.functions.*;

[FE, EF] = S.faceEdges;
edgeFilter = cellfun('length',EF) == 2;
twoFaceEdges = find(edgeFilter);
EE = [EF{edgeFilter}]';
sets = independentSetDecomposition(EE);
score = zeros(1,S.edges.NumObjects);
out = cell(1,length(sets));
better = cell(1,length(sets));
for s=1:length(sets)
    better{s} = twoFaceEdges(sets{s});
end
for s=1:length(sets)
    s
%     score(twoFaceEdges(sets{s})) = deleteEdgesAndScore(S,I, twoFaceEdges(sets{s}));

    out{s} = deleteEdgesAndScore(S,I, better{s});
end
for s=1:length(sets)
    score(twoFaceEdges(sets{s})) = out{s};
end
% NaN?
% score(~edgeFilter) = 0;

end
