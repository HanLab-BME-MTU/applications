function [ A ] = faceAdjacencyFromEdges( S, EF, e )
%faceAdjacencyFromEdges Get face adjacency from edge to face cell map 
% S is a Skeleton class
% EF is a cell array with indices indicated edge number and values indicating adjacent faces
% e are the edges of interest
A = false(S.faces.NumObjects);
EE = [EF{e}]';
for i=1:length(e)
    A(EE(i,:),EE(i,:)) = true;
end
A = A | A';
A = A & ~eye(size(A));

end
