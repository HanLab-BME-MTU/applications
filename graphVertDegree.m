function vDegree = graphVertDegree(edges,nVert)
%GRAPHVERTDEGREE calculates the degree of each vertex in the input skeleton graph
%
% vDegree = skelGraphVertDegree(edges,nVert)
%
%  This function calculates the number of edges connecting to each vertex
%  in the input graph.
%
% Input:
% 
%   edges - The Mx2 list of edges, where each row specifies the starting
%   and ending vertices of that edge
%
%   nVert - The total number of vertices present. Optional. If not input,
%   the maximum value in the edges matrix will be used
%
%
% Output:
%
% vDegree = An nVert x 1 vector specifying the degree of each vertex (the
%           number of edges connecting to it)
%
% Hunter Elliott
% 5/2011
%

if nargin < 1 || isempty(edges) || ndims(edges) ~=2 || ...
        size(edges,2) ~=2 || ~all(isposint(edges(:)))
    error('The first input must be an Mx2 edge matrix containing only positive integers!')
end

if nargin < 2 || isempty(nVert)
    nVert = max(edges(:));
end

vDegree = zeros(nVert,1);
nEdges = size(edges,1);

for j = 1:nEdges    
    if ~any(edges(j,:)==0)
        %Add this edge to the degree count of each vertex it connects
        vDegree(edges(j,:)) = vDegree(edges(j,:)) + 1;
    end
end