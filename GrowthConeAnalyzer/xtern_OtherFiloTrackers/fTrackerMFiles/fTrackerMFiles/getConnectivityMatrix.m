function [connMatrix]=getConnectivityMatrix(vertexesAndNeighbours, labeledImage)

% provides the connectivityMatrix of a Graph
% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


connMatrix=zeros(size(vertexesAndNeighbours, 2));
segmentAreas=regionprops(labeledImage, 'Area');
for itVertexRow=1:size(vertexesAndNeighbours, 2)
    for itVertexCol=itVertexRow+1:size(vertexesAndNeighbours, 2)
        connexions=getConnectivity(vertexesAndNeighbours{itVertexRow},...
            vertexesAndNeighbours{itVertexCol}, labeledImage);
        if max(connexions)~=0
            connMatrix(itVertexRow, itVertexCol)=max(segmentAreas(nonzeros(connexions)).Area);
            connMatrix(itVertexCol, itVertexRow)=connMatrix(itVertexRow, itVertexCol);
        end
    end
end