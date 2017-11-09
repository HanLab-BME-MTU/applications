%[idx1, idx2] = colocalizationLAP(X1, X2, R) determines colocalization between two
% sets of points based on mutual proximity via bipartite graph matching.
%
% Inputs:
%      X1 : first set of points, up to 3D: [x1 y1 z1], NxDim format
%      X2 : second set of points
%
% Outputs:
%    idx1 : index of matched points in the first set 
%    idx2 : index of corresponding points in the second set

% Francois Aguet, 2013

function [idx1, idx2] = colocalizationLAP(X1, X2, R)

D = createSparseDistanceMatrix(X1, X2, R);
[link12, ~] = lap(D, [], [], 1);

n1 = size(X1,1);
n2 = size(X2,1);
link12 = link12(1:n1);
matchIdx = link12<=n2;
idx1 = find(matchIdx);
idx2 = double(link12(matchIdx));
