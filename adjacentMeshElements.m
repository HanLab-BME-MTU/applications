function [iVert,iFace] = adjacentMeshElements(surf,iPts,maxDist)
%ADJACENTMESHELEMENTS finds vertices and faces which are adjacent to the selected vertices on the input surface mesh
% 
%  [iVert,iFace] = adjacentMeshElements(surf,iPts,maxDist)
% 
%  This function finds, for each input vertex, all vertices and faces on
%  the input mesh surface which are within the specified distance of the
%  input vertex, and which are connected to the input vertices by edges
%  which are also within the specified radius of the input points.
%  Therefore, vertices which are within the specified cartesian distance of
%  the input vertices, but are topologically distant along the mesh
%  surface, will not be returned.
% 
%   Input:
% 
%       surf - The structure containing the surface mesh to find adjacent
%       points on, stored in the FV format (faces,vertices) used by patch,
%       isosurface etc.
% 
%       iPts - Mx3 matrix containing the positions of the points to find
%       adjacent mesh vertices for.
%
%       maxDist - A positive scalar or Mx1 vector specifying The maximum
%       distance from the input point to consider adjacent points. If
%       input as a scalar, the same distance is used for all input
%       points, if a vector the same length as iPts, then this specifies a
%       different max distance to use for each selected vertex.
%
%  
%   Output:
%
%       iVert - Mx1 cell array, the m-th element of which contains the
%       indices of the vertices which are within the maxDist of the m-th
%       input vertex.
%
%       iFace - Mx1 cell array, the m-th element of which contains the
%       indices of the faces for which all three vertices are within the
%       maxDist of the m-th input vertex.
%
% Hunter Elliott 
% 4/2011
%

if nargin < 3 || isempty(surf) || isempty(iPts) || isempty(maxDist)        
    error('You must input an FV surface, vertex indices, and a maximum distance!')
end

if ~isfv(surf)
    error('Input surf is not a valid FV surface mesh!')
end

if ~isequal(round(abs(iPts)),iPts)
   error('The input vertex indices must be a positive integer vector or scalar!')
end

nPts = numel(iPts);

if any(maxDist <= 0)
    error('The maxDist input must be a positive scalar or vector!')
elseif numel(maxDist) == 1
    maxDist = ones(nPts,1) .* maxDist;
elseif numel(maxDist) ~= nPts
    error('If maxDist is input as a vector, it must have the same number of elements as the number of input points!')
end

%Make sure we have column vectors
iPts = iPts(:);
maxDist = double(maxDist(:));%Make sure its a double - get strange behavior with other classes


%First, get all the points within specified radii, ignoring mesh topology.
iClose = KDTreeBallQuery(surf.vertices,surf.vertices(iPts,:),maxDist);

iVert  = cell(nPts,1);
if nargout > 1
    iFace = cell(nPts,1);
end

%Now go through each set of points and eliminate those which are not
%connected through mesh edges which are also within the specified radii.
%This prevents topologically-distant points from being returned simply due
%to straight-line proximity.

%Convert the faces to edges to ease the connectivity analysis
edges = vertcat(surf.faces(:,1:2),surf.faces(:,2:3));


for j = 1:nPts
        
    %Get all the edges which connect points close to this input point
    goodEdges = any(bsxfun(@eq,iClose{j}',edges(:,1)),2) | ...
                     any(bsxfun(@eq,iClose{j}',edges(:,2)),2);                     
    
    %Convert these edges into a sparse adjacency matrix (graph)
    maxEdge = max(max(edges(goodEdges,:)));
    edgeGraph = sparse(edges(goodEdges,1),edges(goodEdges,2),...
                1,maxEdge,maxEdge,numel(goodEdges));
            
    %Using only these edges, label the connected components of the
    %resulting graph
    [~,compLabel] = graphconncomp(edgeGraph,'Weak',true);
        
    
    %Return only the vertices which have the same label as this input
    %vertex
    iVert{j} = find(compLabel == compLabel(iPts(j)));
    
    if nargout > 1
        %If requested, return the faces which these vertices compose.
        iFace{j} = find(any(bsxfun(@eq,iVert{j},surf.faces(:,1)),2) & ...
                        any(bsxfun(@eq,iVert{j},surf.faces(:,2)),2) & ...
                        any(bsxfun(@eq,iVert{j},surf.faces(:,3)),2));                                
    end
    
end





