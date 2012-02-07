function edges(obj)

centers = cellfun(@(a) (min(obj.data.points(a(:),:),[],1)+max(obj.data.points(a(:),:),[],1))/2,obj.data.clusters,'UniformOutput',false);

pointsStart = cell2mat(centers(obj.data.edges(:,1),:));
pointsEnd = cell2mat(centers(obj.data.edges(:,2),:));

obj.imaris.displaySegments(pointsStart,pointsEnd,'Display: Edges');

end

