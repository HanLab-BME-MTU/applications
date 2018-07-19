function residuals(obj)

pointsStart = obj.data.points(horzcat(obj.data.clusters{:}),:);
pointsEnd = zeros(obj.data.nClusters,3);

for c=1:obj.data.nClusters
    pointsEnd(c,:) = renderBezier(obj.data.modelBezCP{c},obj.data.modelProj{c});
end

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Residuals');

end



