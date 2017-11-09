function [connMatrix, segmentsIdConnctions]=getConnectivityAndColor(vertexesAndNeighbours, labeledImage)

% provides the connectivityMatrix of a Graph with the 1s replaced by the
% lengh of the shortest connexion segment

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


connMatrix=zeros(size(vertexesAndNeighbours, 2));
segmentsIdConnctions=zeros(size(vertexesAndNeighbours, 2));
segmentAreas=regionprops(labeledImage, 'Area');

for itVertexRow=1:size(vertexesAndNeighbours, 2)
    for itVertexCol=itVertexRow+1:size(vertexesAndNeighbours, 2)
        connexions=getConnectivity(vertexesAndNeighbours{itVertexRow},...
            vertexesAndNeighbours{itVertexCol}, labeledImage);
        if max(connexions)~=0
            noCero=nonzeros(connexions);
            areas=[];
            for itObject=1:size(noCero, 1)
                areas(itObject,1)=noCero(itObject,1);
                areas(itObject,2)=segmentAreas(noCero(itObject,1)).Area;
            end
            [connMatrix(itVertexRow, itVertexCol),minIndex]=min(areas(:,2));
            segmentsIdConnctions(itVertexRow, itVertexCol)=areas(minIndex,1);
            connMatrix(itVertexCol, itVertexRow)=connMatrix(itVertexRow, itVertexCol);
            segmentsIdConnctions(itVertexCol, itVertexRow)=segmentsIdConnctions(itVertexRow, itVertexCol);
        end
    end
end