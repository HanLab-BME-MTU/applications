function sets = independentSetDecomposition(edgeMatrix)
%independentSetDecomposition splits the edgeMatrix into independent sets
% edgeMatrix is a Nx2 matrix indicating an edge from the nodes in the
%   first column to the second column
    sets = cell(1,100);
    E = edgeMatrix;
    % number the edges
    edgeIdx = 1:size(edgeMatrix,1);
    i = 1;
    while(~isempty(E))
        % double the size of the cell array if we expand beyond the size
        if(i > length(sets))
            sets{length(sets)*2} = [];
        end
        M = maxMatching(max(E(:)),E);
        % save the matched edges
        sets{i} = edgeIdx(M);
        % continue matching the unmatched edges
        E = E(~M,:);
        edgeIdx = edgeIdx(~M);
        % track the number of independent sets thus far
        i = i + 1;
    end
    % truncate the sets to only the ones used
    sets = sets(1:i-1);
end
