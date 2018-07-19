function clustersAsChains(obj)

pointsStart = zeros(1,3);
pointsEnd = zeros(1,3);

for c=1:obj.data.nClusters
    points = obj.data.points(obj.data.clusters{c},:);
    pointsStart = [pointsStart;points(1:end-1,:)];
    pointsEnd = [pointsEnd;points(2:end,:)];
end

pointsStart = pointsStart(2:end,:);
pointsEnd = pointsEnd(2:end,:);

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Cluster Chains');

end
