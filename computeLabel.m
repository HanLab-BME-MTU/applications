function labels = computeLabel(E, N)
% label the vertices of a graph defined by edges E. N is the number of
% vertices in the graph. labels(v) == connected component label in which v
% belongs to. labels(v) is in [1...nCC] where nCC is the number of
% connected components in the graph.

G = sparse(E(:,1),E(:,2),true(size(E,1),1), N, N, size(E,1));

CC = cell(N,1);

marked = false(N,1);

for iTrack = 1:N
    if ~marked(iTrack)
        [CC{iTrack} marked] = rec(G, iTrack, marked);
    end
end

ind = cellfun(@isempty, CC);
CC = CC(~ind);

labels = zeros(N,1);

for iLabel = 1:numel(CC)
    labels(ismember(1:N,CC{iLabel})) = iLabel;
end

function [inCC marked] = rec(G, i, marked)
marked(i) = true;

jdx = find(G(i,:));

nInCC = numel(jdx);
inCC= cell(1,nInCC+1);
inCC{1} = i;

for jj = 1:numel(jdx)
    j = jdx(jj);
    
    if (~marked(j))
        [inCC{jj+1} marked] = rec(G, j, marked);
    end
end

ind = cellfun(@isempty, inCC);
inCC = horzcat(inCC{~ind});
