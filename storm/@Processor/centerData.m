function centerData(obj)
% Compute the center of the point cloud
center = (max(obj.data.points,[],1)+min(obj.data.points,[],1))/2;

% Translate all the points
obj.data.points = obj.data.points+repmat(-center,size(obj.data.points,1),1);
end